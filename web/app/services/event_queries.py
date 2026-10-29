"""Read-only, paginated access to a completed analysis event index."""
from __future__ import annotations

from contextlib import contextmanager
import json
from pathlib import Path
import re
import sqlite3

from reaxnet.parser import canonical_symbol


class EventIndexUnavailable(ValueError):
    """A completed task has no readable index compatible with this UI."""


def normalize_formula(value: str) -> str:
    """Accept plain formulas only, returning an element-order-independent key."""
    value = value.strip()
    tokens = list(re.finditer(r"([A-Z][a-z]?)([0-9]*)", value))
    if not value or ''.join(match.group() for match in tokens) != value:
        raise ValueError('分子式须为元素符号和正整数，例如 C2H6O')
    counts = {}
    for token in tokens:
        symbol = canonical_symbol(token.group(1))
        count = int(token.group(2) or '1')
        if count <= 0:
            raise ValueError('分子式中的原子数必须为正整数')
        counts[symbol] = counts.get(symbol, 0) + count
    order = (['C'] + (['H'] if 'H' in counts else []) + sorted(counts.keys() - {'C', 'H'})
             if 'C' in counts else sorted(counts))
    return ''.join(symbol + (str(counts[symbol]) if counts[symbol] != 1 else '') for symbol in order)


@contextmanager
def _open_index(path: Path):
    if not path.is_file():
        raise EventIndexUnavailable('此任务没有可用的事件索引，请使用当前版本重新分析。')
    db = None
    try:
        db = sqlite3.connect(path.resolve().as_uri() + '?mode=ro', uri=True, timeout=2)
        db.row_factory = sqlite3.Row
        db.execute('PRAGMA query_only = ON')
        yield db
    except sqlite3.Error as exc:
        raise EventIndexUnavailable('事件索引无法读取或版本不兼容，请重新分析此任务。') from exc
    finally:
        if db is not None:
            db.close()


def list_events(path: Path, *, offset=0, limit=50, q='', formula='', species_key='',
                start_step=None, end_step=None, status='all') -> dict:
    if start_step is not None and end_step is not None and start_step > end_step:
        raise ValueError('起始步数不能大于结束步数')
    if status not in {'all', 'accepted', 'raw', 'recrossing'}:
        raise ValueError('未知事件状态')
    where, params = [], []
    if q:
        where.append('''(instr(lower(e.event_id), lower(?)) > 0
            OR instr(lower(e.reaction), lower(?)) > 0
            OR EXISTS (SELECT 1 FROM event_species s WHERE s.event_id=e.event_id
                AND (instr(lower(s.formula), lower(?)) > 0
                     OR instr(lower(s.species_key), lower(?)) > 0)))''')
        params.extend([q] * 4)
    if formula:
        where.append('EXISTS (SELECT 1 FROM event_species s WHERE s.event_id=e.event_id AND s.formula=?)')
        params.append(normalize_formula(formula))
    if species_key:
        where.append('EXISTS (SELECT 1 FROM event_species s WHERE s.event_id=e.event_id AND s.species_key=?)')
        params.append(species_key)
    # A saved event belongs to an interval; range filters include overlaps.
    if start_step is not None:
        where.append('e.timestep >= ?')
        params.append(start_step)
    if end_step is not None:
        where.append('e.previous_timestep <= ?')
        params.append(end_step)
    if status != 'all':
        where.append('e.accepted = ?')
        params.append(0 if status == 'recrossing' else 1)
    clause = ' WHERE ' + ' AND '.join(where) if where else ''
    with _open_index(path) as db:
        total = db.execute('SELECT COUNT(*) FROM events e' + clause, params).fetchone()[0]
        cursor = db.execute('''SELECT e.event_id, e.previous_timestep, e.timestep,
            e.time_start_ps, e.time_end_ps, e.reaction, e.accepted, e.filter_reason
            FROM events e''' + clause +
            ' ORDER BY e.timestep, e.event_id LIMIT ? OFFSET ?', [*params, limit, offset])
        rows = []
        for row in cursor:
            result = dict(row)
            result['accepted'] = bool(result['accepted'])
            rows.append(result)
    return {'rows': rows, 'total': total, 'offset': offset, 'limit': limit,
            'has_more': offset + len(rows) < total}


def get_event(path: Path, event_id: str) -> dict | None:
    with _open_index(path) as db:
        record = db.execute('SELECT payload FROM events WHERE event_id=?', (event_id,)).fetchone()
        if record is None:
            return None
        try:
            payload = json.loads(record['payload'])
            if not isinstance(payload, dict):
                raise ValueError('event payload is not an object')
        except (ValueError, TypeError) as exc:
            raise EventIndexUnavailable('事件索引中的记录损坏，请重新分析。') from exc
        return payload
