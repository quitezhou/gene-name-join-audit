"""Core logic: AnnData-style name assignment, two-way join audit, and ID recovery.

Standard library only. The algorithms are those used in the accompanying
manuscript (``audit_join.py`` and ``recover_ids.py``), reorganised as a library.
"""

from __future__ import annotations

import re
from bisect import bisect_left, bisect_right
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from typing import Iterable, Iterator, Sequence


# ---------------------------------------------------------------------------
# Name assignment
# ---------------------------------------------------------------------------

def assign_keys(symbols: Sequence[str], separator: str = "-") -> Iterator[str]:
    """Yield unique keys the way ``AnnData.var_names_make_unique`` does.

    The first occurrence of a name keeps it; later occurrences receive
    ``<name><separator><n>`` with the smallest n not already used and not
    present among the original names.
    """
    reserved = set(symbols)
    seen: set[str] = set()
    last: Counter[str] = Counter()
    for name in symbols:
        if name not in seen:
            seen.add(name)
            yield name
            continue
        while True:
            last[name] += 1
            candidate = f"{name}{separator}{last[name]}"
            if candidate not in reserved:
                reserved.add(candidate)
                yield candidate
                break


# ---------------------------------------------------------------------------
# Gene axes
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Row:
    id: str
    symbol: str
    key: str


class Axis:
    """One side of a join: stable ID, original symbol and final key per row.

    ``entries`` is a sequence of ``(id, symbol, key)``. If every key is
    ``None`` the keys are generated with :func:`assign_keys` in the given order.
    """

    def __init__(self, entries: Sequence[tuple[str, str, str | None]], separator: str = "-",
                 *, symbols_known: bool = True):
        if not entries:
            raise ValueError("No gene records supplied")
        keys: Iterable[str | None] = [r[2] for r in entries]
        if all(k is None for k in keys):
            keys = assign_keys([r[1] for r in entries], separator)
        self.by_key: dict[str, Row] = {}
        self.by_id: dict[str, list[Row]] = defaultdict(list)
        self.by_symbol: dict[str, set[str]] = defaultdict(set)
        self.multiplicity: Counter[str] = Counter()
        self.symbols_known = symbols_known
        for (fid, symbol, _), key in zip(entries, keys):
            if not all(isinstance(v, str) and v.strip() for v in (fid, symbol, key)) or key in self.by_key:
                raise ValueError(f"Missing ID, symbol or key, or duplicate final key: {(fid, symbol, key)}")
            row = Row(fid, symbol, key)
            self.by_key[key] = row
            self.by_id[fid].append(row)
            self.by_symbol[symbol].add(fid)
            self.multiplicity[symbol] += 1

    def __len__(self) -> int:
        return len(self.by_key)

    def use_source_groups(self, source: "Axis") -> None:
        """Keep this axis's rows and keys, but judge duplicate families on the
        complete source table it was filtered from."""
        if self.symbols_known:
            for row in self.by_key.values():
                if row.id not in source.by_symbol.get(row.symbol, ()):
                    raise ValueError(f"Analysed feature absent from supplied source table: {row.id}, {row.symbol}")
        else:
            entries = []
            for row in self.by_key.values():
                symbols = {s.symbol for s in source.by_id.get(row.id, ())}
                if len(symbols) != 1:
                    raise ValueError(f"Source table must provide exactly one original symbol for {row.id}")
                entries.append((row.id, symbols.pop(), row.key))
            restored = Axis(entries)
            self.by_key, self.by_id = restored.by_key, restored.by_id
            self.symbols_known = True
        self.by_symbol = source.by_symbol
        self.multiplicity = source.multiplicity


# ---------------------------------------------------------------------------
# Join audit
# ---------------------------------------------------------------------------

EVENT_FIELDS = ("event", "class", "left_id", "right_id", "left_symbol", "right_symbol",
                "left_key", "right_key", "left_multiplicity", "right_multiplicity",
                "cross_table_family_id_overlap", "left_key_absent_right", "right_key_absent_left")


@dataclass
class Report:
    """Result of :func:`audit`.

    ``same_name_different_id`` and ``same_id_different_name`` are the two
    checks recommended before joining. ``ok`` is True only when both are zero.
    """

    summary: dict
    events: list[dict] = field(default_factory=list)

    @property
    def same_name_different_id(self) -> int:
        return (self.summary["different_ids_duplicate"] + self.summary["different_ids_singleton"]
                + self.summary.get("different_ids_unclassified", 0))

    @property
    def same_id_different_name(self) -> int:
        return self.summary["shared_id_different_keys"]

    @property
    def ok(self) -> bool:
        return self.same_name_different_id == 0 and self.same_id_different_name == 0

    def to_dict(self) -> dict:
        return {"ok": self.ok, "same_name_different_id": self.same_name_different_id,
                "same_id_different_name": self.same_id_different_name, "summary": self.summary}

    def __repr__(self) -> str:
        s = self.summary
        return (f"Report(ok={self.ok}, shared_names={s['shared_keys']}, shared_ids={s['shared_ids']}, "
                f"same_name_different_id={self.same_name_different_id}, "
                f"same_id_different_name={self.same_id_different_name})")


def audit(left: Axis, right: Axis) -> Report:
    """Compare every shared key with its IDs, and every shared ID with its keys."""
    shared_keys = left.by_key.keys() & right.by_key.keys()
    shared_ids = left.by_id.keys() & right.by_id.keys()
    summary = dict(shared_keys=len(shared_keys), shared_ids=len(shared_ids),
                   different_ids_duplicate=0, different_ids_singleton=0, duplicate_with_family_id_overlap=0,
                   shared_ids_nonunique_in_either_axis=0, shared_id_different_keys=0,
                   shared_id_changed_symbol_key_preserved=0,
                   singly_occurring_renamed_ids=0, singleton_renames_both_keys_absent_opposite=0)
    events: list[dict] = []
    both_symbols_known = left.symbols_known and right.symbols_known

    def add(event, kind, a, b, overlap):
        events.append(dict(zip(EVENT_FIELDS, (
            event, kind, a.id, b.id, a.symbol if left.symbols_known else "",
            b.symbol if right.symbols_known else "", a.key, b.key,
            left.multiplicity[a.symbol] if left.symbols_known else "",
            right.multiplicity[b.symbol] if right.symbols_known else "", overlap,
            a.key not in right.by_key, b.key not in left.by_key))))

    for key in sorted(shared_keys):
        a, b = left.by_key[key], right.by_key[key]
        if a.id == b.id:
            continue
        duplicate = left.multiplicity[a.symbol] > 1 or right.multiplicity[b.symbol] > 1
        if not both_symbols_known and not duplicate:
            summary["different_ids_unclassified"] = summary.get("different_ids_unclassified", 0) + 1
            add("same_key_different_id", "unknown_symbol_provenance", a, b, "")
            continue
        overlap = both_symbols_known and (a.id in right.by_symbol.get(a.symbol, ())
                                          or b.id in left.by_symbol.get(b.symbol, ()))
        summary["different_ids_duplicate" if duplicate else "different_ids_singleton"] += 1
        summary["duplicate_with_family_id_overlap"] += int(duplicate and overlap)
        add("same_key_different_id", "duplicate_associated" if duplicate else "singleton_on_both",
            a, b, overlap if both_symbols_known else "")
    for fid in sorted(shared_ids):
        if len(left.by_id[fid]) != 1 or len(right.by_id[fid]) != 1:
            summary["shared_ids_nonunique_in_either_axis"] += 1
            continue
        a, b = left.by_id[fid][0], right.by_id[fid][0]
        if a.key == b.key:
            if both_symbols_known and a.symbol != b.symbol:
                summary["shared_id_changed_symbol_key_preserved"] += 1
                add("same_id_changed_symbol_key_preserved", "annotation_change_only", a, b, "")
            continue
        summary["shared_id_different_keys"] += 1
        unique_rename = (both_symbols_known and left.multiplicity[a.symbol] == right.multiplicity[b.symbol] == 1
                         and a.symbol != b.symbol)
        summary["singly_occurring_renamed_ids"] += int(unique_rename)
        summary["singleton_renames_both_keys_absent_opposite"] += int(
            unique_rename and a.key not in right.by_key and b.key not in left.by_key)
        kind = "other_key_change" if both_symbols_known else "unknown_symbol_provenance"
        add("same_id_different_keys", "singleton_rename" if unique_rename else kind, a, b, "")
    return Report(summary, events)


# ---------------------------------------------------------------------------
# Recovery when IDs were lost
# ---------------------------------------------------------------------------

_SUFFIX = re.compile(r"(.+)-[1-9][0-9]*")


class Reference:
    """Complete source feature table used to bound the identity of final keys."""

    def __init__(self, records: Sequence[tuple[str, str]]):
        self.ids: list[str] = []
        self.positions: dict[str, list[int]] = defaultdict(list)
        for row, (identifier, symbol) in enumerate(records):
            if not all(isinstance(v, str) and v.strip() for v in (identifier, symbol)):
                raise ValueError("Empty source identifier or symbol")
            self.ids.append(identifier)
            self.positions[symbol].append(row)
        if not self.ids or len(set(self.ids)) != len(self.ids):
            raise ValueError("Source must contain nonempty, one-to-one identifiers")

    def possible_rows(self, label: str) -> list[int]:
        # Keep real source names that happen to end in e.g. "-1". A generated
        # suffix requires at least two identical original names on the source.
        possible = self.positions.get(label, [])
        suffix = _SUFFIX.fullmatch(label)
        if suffix:
            base = self.positions.get(suffix[1], [])
            if len(base) > 1:
                return sorted([*possible, *base])
        return possible

    def candidates(self, labels: Sequence[str], preserve_order: bool = False) -> dict[str, tuple[str, ...]]:
        """Return, for each final key, the source IDs it could denote.

        With ``preserve_order=True`` (use only when provenance shows that rows
        kept their source order) one forward and one backward scan remove
        positions that no order-preserving placement allows. A single
        candidate identifies the record under these assumptions.
        """
        if (not labels or not all(isinstance(key, str) and key.strip() for key in labels)
                or len(set(labels)) != len(labels)):
            raise ValueError("Supply the complete ordered list of nonempty, unique final keys")
        rows = [self.possible_rows(key) for key in labels]
        for key, choices in zip(labels, rows):
            if not choices:
                raise ValueError(f"No compatible source symbol for {key!r}")
        if preserve_order:
            earliest = []
            previous = -1
            for key, choices in zip(labels, rows):
                index = bisect_right(choices, previous)
                if index == len(choices):
                    raise ValueError(f"Order check failed at {key!r}: these keys cannot come from this source "
                                     "in this order. Do not use --order-preserved.")
                previous = choices[index]
                earliest.append(previous)
            following = len(self.ids)
            for j in range(len(rows) - 1, -1, -1):
                choices = rows[j]
                end = bisect_left(choices, following)
                start = bisect_right(choices, earliest[j - 1] if j else -1)
                rows[j] = choices[start:end]
                following = choices[end - 1]
        return {label: tuple(self.ids[i] for i in choices) for label, choices in zip(labels, rows)}


def recover(source: Sequence[tuple[str, str]], keys: Sequence[str],
            preserve_order: bool = False) -> dict[str, tuple[str, ...]]:
    """Convenience wrapper around :class:`Reference`."""
    return Reference(source).candidates(keys, preserve_order)
