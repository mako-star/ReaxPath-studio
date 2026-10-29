"""Disk-backed, finalized candidate-event evidence for pagination and inspection."""
from __future__ import annotations

import json
import sqlite3
from pathlib import Path


class EventIndex:
    def __init__(self, path):
        self.path = Path(path)
        self.db = None

    def __enter__(self):
        self.db = sqlite3.connect(self.path)
        self.db.executescript("""
            CREATE TABLE IF NOT EXISTS events (
                event_id TEXT PRIMARY KEY, previous_timestep INTEGER, timestep INTEGER,
                time_start_ps REAL, time_end_ps REAL, reaction TEXT,
                accepted INTEGER, filter_reason TEXT, payload TEXT);
            CREATE INDEX IF NOT EXISTS events_by_step ON events(timestep,event_id);
            CREATE INDEX IF NOT EXISTS events_by_accepted ON events(accepted,timestep,event_id);
            CREATE TABLE IF NOT EXISTS event_species (event_id TEXT, species_key TEXT, formula TEXT,
                PRIMARY KEY (event_id,species_key));
            CREATE INDEX IF NOT EXISTS event_species_by_key ON event_species(species_key,event_id);
            CREATE INDEX IF NOT EXISTS event_species_by_formula ON event_species(formula,event_id);
            DELETE FROM event_species;
            DELETE FROM events;
        """)
        return self

    def write(self, event, registry):
        from reaxnet.outputs import event_to_payload
        payload = event_to_payload(event, registry)
        self.db.execute('INSERT INTO events VALUES (?,?,?,?,?,?,?,?,?)', (
            event.event_id,event.previous_timestep,event.timestep,
            payload['time_start_ps'],payload['time_end_ps'],payload['reaction'],
            int(event.accepted),event.filter_reason,
            json.dumps(payload,ensure_ascii=False,separators=(',',':'),allow_nan=False)))
        molecules = payload['reactants'] + payload['products']
        self.db.executemany('INSERT OR IGNORE INTO event_species VALUES (?,?,?)',
            [(event.event_id,m['species_key'],m['formula']) for m in molecules])

    def __exit__(self, exc_type, exc, traceback):
        try:
            if exc_type is None: self.db.commit()
            else: self.db.rollback()
        finally:
            self.db.close()
