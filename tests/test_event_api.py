from __future__ import annotations

import csv
from contextlib import closing
import json
from pathlib import Path
import sqlite3
import sys
import tempfile
import unittest
from unittest.mock import patch

import httpx

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / 'web') not in sys.path:
    sys.path.insert(0, str(ROOT / 'web'))

from app.main import create_app
from app.services import analyzer, job_manager, storage


class EventApiTests(unittest.IsolatedAsyncioTestCase):
    async def asyncSetUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.patches = [patch.object(job_manager, 'JOB_ROOT', self.root / 'jobs'),
                        patch.object(analyzer, 'JOB_ROOT', self.root / 'jobs'),
                        patch.object(storage, 'UPLOAD_ROOT', self.root / 'uploads')]
        for active in self.patches:
            active.start()
        job_manager._singleton = None
        app = create_app()
        self.client = httpx.AsyncClient(transport=httpx.ASGITransport(app=app,
                                      client=('127.0.0.1', 43210)), base_url='http://127.0.0.1')
        self.manager = job_manager.get_job_manager()
        self.index = self.root / 'event_index.sqlite'
        self._make_index(self.index)
        await self._job('events-test-001', {'event_index': str(self.index)})

    async def asyncTearDown(self):
        await self.client.aclose()
        job_manager._singleton = None
        for active in reversed(self.patches):
            active.stop()
        self.temporary.cleanup()

    async def _job(self, job_id, files):
        await self.manager.create(job_id, 'unused', 'sample.bonds', {}, {})
        await self.manager.update(job_id, status=job_manager.JobStatus.COMPLETED, result_files=files)

    def _make_index(self, path, count=3):
        with closing(sqlite3.connect(path)) as db, db:
            db.executescript('''
                CREATE TABLE events (event_id TEXT PRIMARY KEY, previous_timestep INTEGER,
                    timestep INTEGER, time_start_ps REAL, time_end_ps REAL, reaction TEXT,
                    accepted INTEGER, filter_reason TEXT, payload TEXT);
                CREATE TABLE event_species (event_id TEXT, species_key TEXT, formula TEXT);
                CREATE INDEX events_step ON events(timestep);
                CREATE INDEX event_formula ON event_species(formula, event_id);
            ''')
            for i in range(count):
                event_id = f'E{i + 1:06d}'
                formula = ('CH4', 'CH3', 'C2H6O')[i]
                molecule = {'instance_id': f'F{i:08d}:M00000000', 'species_key': f'S{i}',
                            'structure_id': f'{formula}-001', 'formula': formula,
                            'atoms': [{'atom_id': 1, 'element': 'C'}], 'connections': [],
                            'smiles': '', 'inchi': '', 'inchikey': '',
                            'identity_status': 'disabled', 'message': ''}
                payload = {'event_id': event_id, 'previous_timestep': i * 10,
                           'timestep': (i + 1) * 10, 'time_start_ps': i / 10,
                           'time_end_ps': (i + 1) / 10, 'reaction': f'{formula} → X',
                           'accepted': i != 1, 'filter_reason': 'return' if i == 1 else '',
                           'reactants': [molecule], 'products': [molecule],
                           'event_type': 'rearrangement', 'broken_bonds': [[1, 2]],
                           'formed_bonds': [[1, 3]], 'changed_atom_ids': [1, 2, 3],
                           'conservation': {'atom_ids_equal': True}}
                db.execute('INSERT INTO events VALUES (?,?,?,?,?,?,?,?,?)',
                           (event_id, i * 10, (i + 1) * 10, i / 10, (i + 1) / 10,
                            payload['reaction'], int(i != 1), payload['filter_reason'], json.dumps(payload)))
                # A formula appearing on both sides must not duplicate an event.
                db.executemany('INSERT INTO event_species VALUES (?,?,?)',
                               [(event_id, f'S{i}', formula), (event_id, f'S{i}', formula)])

    async def test_exact_pagination_and_summary_without_loading_all_details(self):
        first = await self.client.get('/api/jobs/events-test-001/events', params={'limit': 2})
        self.assertEqual(first.status_code, 200, first.text)
        data = first.json()
        self.assertEqual(data['total'], 3)
        self.assertTrue(data['has_more'])
        self.assertEqual([r['event_id'] for r in data['rows']], ['E000001', 'E000002'])
        self.assertNotIn('reactants', data['rows'][0])
        last = (await self.client.get('/api/jobs/events-test-001/events',
                                      params={'offset': 2, 'limit': 1})).json()
        self.assertEqual(last['rows'][0]['event_id'], 'E000003')
        self.assertFalse(last['has_more'])
        beyond = (await self.client.get('/api/jobs/events-test-001/events',
                                        params={'offset': 50})).json()
        self.assertEqual(beyond['rows'], [])
        self.assertEqual(beyond['total'], 3)

    async def test_combined_filters_formula_normalization_and_literal_search(self):
        response = await self.client.get('/api/jobs/events-test-001/events',
            params={'formula': 'H6C2O', 'species_key': 'S2', 'q': 'c2h6o', 'status': 'accepted',
                    'start_step': 25, 'end_step': 25})
        self.assertEqual(response.status_code, 200, response.text)
        self.assertEqual([r['event_id'] for r in response.json()['rows']], ['E000003'])
        for query in ("' OR 1=1 --", '%', '_'):
            response = await self.client.get('/api/jobs/events-test-001/events', params={'q': query})
            self.assertEqual(response.json()['total'], 0)
        for state, count in (('all', 3), ('accepted', 2), ('raw', 2), ('recrossing', 1)):
            response = await self.client.get('/api/jobs/events-test-001/events', params={'status': state})
            self.assertEqual(response.json()['total'], count)

    async def test_event_list_does_not_load_or_decode_full_evidence(self):
        # Invalid evidence JSON must not affect the scalar-only listing. Detail
        # remains responsible for checking and decoding the full record.
        with closing(sqlite3.connect(self.index)) as db, db:
            db.execute("UPDATE events SET payload=? WHERE event_id='E000001'", ('{' * 1000000,))
        listing = await self.client.get('/api/jobs/events-test-001/events', params={'limit': 1})
        self.assertEqual(listing.status_code, 200, listing.text)
        self.assertEqual(listing.json()['rows'][0]['event_id'], 'E000001')
        self.assertNotIn('payload', listing.json()['rows'][0])
        detail = await self.client.get('/api/jobs/events-test-001/events/E000001')
        self.assertEqual(detail.status_code, 409)

    async def test_detail_and_single_event_export_preserve_evidence(self):
        path = '/api/jobs/events-test-001/events/E000001'
        detail = await self.client.get(path)
        self.assertEqual(detail.status_code, 200, detail.text)
        self.assertEqual(detail.json()['reactants'][0]['atoms'], [{'atom_id': 1, 'element': 'C'}])
        self.assertEqual(detail.json()['broken_bonds'], [[1, 2]])
        exported = await self.client.get(path + '/export')
        self.assertEqual(exported.json(), detail.json())
        self.assertIn('attachment', exported.headers['content-disposition'])
        self.assertEqual(exported.headers['content-type'], 'application/json')
        missing = await self.client.get('/api/jobs/events-test-001/events/E999999')
        self.assertEqual(missing.status_code, 404)

    async def test_empty_old_missing_and_incompatible_indices(self):
        empty = self.root / 'empty.sqlite'
        self._make_index(empty, 0)
        await self._job('events-empty-001', {'event_index': str(empty)})
        response = await self.client.get('/api/jobs/events-empty-001/events')
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.json()['total'], 0)
        self.assertEqual(response.json()['rows'], [])
        await self._job('events-legacy-01', {})
        response = await self.client.get('/api/jobs/events-legacy-01/events')
        self.assertEqual(response.status_code, 409)
        self.assertIn('重新分析', response.json()['detail'])
        await self._job('events-missing-1', {'event_index': str(self.root / 'gone.sqlite')})
        self.assertEqual((await self.client.get('/api/jobs/events-missing-1/events')).status_code, 409)
        incompatible = self.root / 'incompatible.sqlite'
        with closing(sqlite3.connect(incompatible)) as db, db:
            db.execute('CREATE TABLE old_events (id TEXT)')
        await self._job('events-old-db-01', {'event_index': str(incompatible)})
        self.assertEqual((await self.client.get('/api/jobs/events-old-db-01/events')).status_code, 409)
        self.assertEqual((await self.client.get('/api/jobs/no-such-job/events')).status_code, 404)

    async def test_invalid_ranges_and_parameters_are_rejected(self):
        for params in ({'start_step': 30, 'end_step': 10}, {'limit': 0}, {'limit': 1001},
                       {'offset': -1}, {'start_step': 'nan'}, {'end_step': 'inf'},
                       {'formula': 'C0H2'}, {'formula': 'Qq2'}, {'formula': 'C2H6O;DROP'},
                       {'status': 'anything'}):
            with self.subTest(params=params):
                response = await self.client.get('/api/jobs/events-test-001/events', params=params)
                self.assertEqual(response.status_code, 422, response.text)
        for params in ({'timestep_scale': 'nan'}, {'timestep_scale': 'inf'},
                       {'recrossing_window': '-inf'}, {'expected_interval': 0}):
            response = await self.client.post('/api/jobs/local', params={
                'bonds_path': str(ROOT / 'examples/ch4_dissociation.bonds'),
                'element_mapping': '1:C,2:H', **params})
            self.assertEqual(response.status_code, 422, response.text)

    async def test_csv_preview_offset_has_more_and_download_mimetypes(self):
        csv_path = self.root / 'registry.csv'
        with csv_path.open('w', newline='', encoding='utf-8') as handle:
            writer = csv.DictWriter(handle, fieldnames=['species_key'])
            writer.writeheader()
            writer.writerows([{'species_key': f'S{i}'} for i in range(3)])
        jsonl = self.root / 'details.jsonl'
        jsonl.write_text('{}\n', encoding='utf-8')
        await self.manager.update('events-test-001', result_files={
            'molecule_registry': str(csv_path), 'event_index': str(self.index), 'event_details': str(jsonl)})
        response = await self.client.get('/api/jobs/events-test-001/preview/molecule_registry',
                                         params={'offset': 1, 'limit': 2})
        self.assertEqual(response.json()['rows'], [{'species_key': 'S1'}, {'species_key': 'S2'}])
        self.assertFalse(response.json()['has_more'])
        self.assertFalse(response.json()['possibly_truncated'])
        first = await self.client.get('/api/jobs/events-test-001/preview/molecule_registry', params={'limit': 2})
        self.assertTrue(first.json()['has_more'])
        for name, expected in [('event_index', 'application/vnd.sqlite3'),
                                ('event_details', 'application/x-ndjson')]:
            download = await self.client.get(f'/api/jobs/events-test-001/outputs/{name}')
            self.assertEqual(download.status_code, 200)
            self.assertEqual(download.headers['content-type'], expected)

    async def test_real_toy_job_exposes_indexed_event_and_conserved_atoms(self):
        response = await self.client.post('/api/jobs/local', params={
            'bonds_path': str(ROOT / 'examples/ch4_dissociation.bonds'),
            'element_mapping': '1:C,2:H', 'timestep_scale': 0.001,
            'generate_identifiers': 'false', 'expected_interval': 10,
        })
        self.assertEqual(response.status_code, 201, response.text)
        job_id = response.json()['job_id']
        job = (await self.client.get(f'/api/jobs/{job_id}')).json()
        self.assertEqual(job['status'], 'completed', job)
        self.assertEqual(job['config']['expected_interval'], 10)
        self.assertIn('event_index', job['result_files'])
        listing = await self.client.get(f'/api/jobs/{job_id}/events', params={'formula': 'CH4'})
        self.assertEqual(listing.status_code, 200, listing.text)
        self.assertEqual(listing.json()['total'], 1)
        event_id = listing.json()['rows'][0]['event_id']
        detail = await self.client.get(f'/api/jobs/{job_id}/events/{event_id}')
        self.assertEqual(detail.status_code, 200, detail.text)
        event = detail.json()
        self.assertEqual(event['previous_timestep'], 0)
        self.assertEqual(event['timestep'], 10)
        self.assertEqual(event['time_start_ps'], 0)
        self.assertEqual(event['time_end_ps'], 0.01)
        self.assertEqual(event['conservation'], {'atom_ids': True, 'elements': True})
        left = {a['atom_id'] for m in event['reactants'] for a in m['atoms']}
        right = {a['atom_id'] for m in event['products'] for a in m['atoms']}
        self.assertEqual(left, {1, 2, 3, 4, 5})
        self.assertEqual(left, right)
        self.assertEqual(event['broken_bonds'], [[1, 5]])
        exported = await self.client.get(f'/api/jobs/{job_id}/events/{event_id}/export')
        self.assertEqual(exported.json(), event)


if __name__ == '__main__':
    unittest.main()
