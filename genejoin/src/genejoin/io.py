"""Readers for feature tables, key traces, 10x outputs and AnnData objects."""

from __future__ import annotations

import csv
import gzip
import re
from pathlib import Path
from typing import Sequence

from .core import Axis, assign_keys

GENE_TYPE = "Gene Expression"
ID_COLUMNS = ("gene_ids", "gene_id", "ensembl_id", "ensembl_ids", "gene_ensembl_id", "feature_id", "id", "ID")
SYMBOL_COLUMNS = ("gene_symbols", "gene_symbol", "feature_name", "gene_name", "symbol", "name")
ENSEMBL_GENE_ID = re.compile(r"ENS[A-Z]*G[0-9]+(?:\.[0-9]+)?")


def _open(path: Path):
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    return opener(path, "rt", encoding="utf-8-sig", newline="")


def read_feature_table(path, mixed_modalities: str = "refuse") -> list[tuple[str, str]]:
    """Read a headerless 10x-style table: ID, symbol [, feature type].

    ``mixed_modalities``:
      * ``"refuse"`` (default): error if non-gene rows are present, because the
        workflow may have named all features together.
      * ``"genes-only"``: drop non-gene rows before anything else.
    """
    if mixed_modalities not in ("refuse", "genes-only"):
        raise ValueError("mixed_modalities must be 'refuse' or 'genes-only'")
    path = resolve_input(path)
    rows, other = [], 0
    with _open(path) as handle:
        for number, fields in enumerate(csv.reader(handle, delimiter="\t"), 1):
            if not fields:
                continue
            if len(fields) not in (2, 3):
                raise ValueError(f"{path}: line {number} has {len(fields)} columns; expected ID, symbol [, type]")
            if len(fields) == 3 and fields[2] != GENE_TYPE:
                other += 1
                continue
            rows.append((fields[0], fields[1]))
    if other and mixed_modalities == "refuse":
        raise ValueError(
            f"{path}: {other} non-gene rows found. If your workflow named all features together, export the "
            "actual final keys and use a keyed table; otherwise pass --mixed-modalities genes-only.")
    return rows


def read_keyed_table(path) -> list[tuple[str, str, str]]:
    """Read a TSV with header ``stable_id, symbol, final_key`` (the actual keys used)."""
    with _open(path) as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, None)
        if header != ["stable_id", "symbol", "final_key"]:
            raise ValueError(f"{path}: expected header stable_id<TAB>symbol<TAB>final_key")
        out = []
        for number, fields in enumerate(reader, 2):
            if len(fields) != 3:
                raise ValueError(f"{path}: line {number} has {len(fields)} columns; expected 3")
            out.append(tuple(fields))
        return out


def read_keys(path) -> list[str]:
    """One final key per line, no header."""
    with _open(path) as handle:
        return [line.rstrip("\r\n") for line in handle]


def _pick(columns, preferred, explicit):
    if explicit:
        if explicit not in columns:
            raise ValueError(f"Column {explicit!r} not found in .var (available: {list(columns)})")
        return explicit
    for name in preferred:
        if name in columns:
            return name
    return None


def _text_values(values, label):
    result = list(values)
    if any(not isinstance(value, str) or not value.strip() for value in result):
        raise ValueError(f"{label} must contain a nonempty text value for every feature")
    return result


def _gene_mask(var):
    """Read CELLxGENE gene/spike-in labels and 10x assay modalities."""
    mask = [True] * len(var)
    modality_columns = ("feature_types", "feature_type")
    if "feature_biotype" in var.columns:
        biotypes = list(var["feature_biotype"])
        cellxgene_labels = {"gene", "spike-in"}
        if any(isinstance(value, str) and value in cellxgene_labels for value in biotypes):
            if not all(isinstance(value, str) and value in cellxgene_labels for value in biotypes):
                raise ValueError("CELLxGENE feature_biotype must contain gene or spike-in for every row")
            mask = [value == "gene" for value in biotypes]
            # CELLxGENE's singular feature_type is a gene classification such as protein_coding.
            modality_columns = ("feature_types",)
    column = _pick(var.columns, modality_columns, None)
    if column is not None:
        modalities = _text_values(var[column], column)
        mask = [keep and modality == GENE_TYPE for keep, modality in zip(mask, modalities)]
    return mask


def axis_from_anndata(adata, id_column: str | None = None, symbol_column: str | None = None,
                      mixed_modalities: str = "refuse") -> Axis:
    """Build an :class:`Axis` from an AnnData-like object without reading the matrix.

    Uses ``adata.var_names`` as the final keys, an Ensembl ID column from
    ``adata.var`` and, if present, an original-symbol column. Without a symbol
    column the two identity checks still run, but mechanism attribution may
    require a complete source table.
    """
    var = adata.var
    if mixed_modalities not in ("refuse", "genes-only"):
        raise ValueError("mixed_modalities must be 'refuse' or 'genes-only'")
    keys = list(adata.var_names)
    mask = _gene_mask(var)
    if not all(mask):
        if mixed_modalities == "refuse":
            raise ValueError("Non-gene features present. Use --mixed-modalities genes-only to audit gene rows.")
        if not any(mask):
            raise ValueError("No gene records remain after filtering non-gene features")
        var = var.loc[mask]
        keys = [key for key, keep in zip(keys, mask) if keep]
    keys = _text_values(keys, ".var_names")
    idc = _pick(var.columns, ID_COLUMNS, id_column)
    if idc is None:
        if keys and all(ENSEMBL_GENE_ID.fullmatch(k) for k in keys):
            ids = keys
        else:
            raise ValueError(
                "No stable ID column found in .var. Pass id_column=..., or, if the IDs were lost, use "
                "`genejoin recover` with the original feature table.")
    else:
        ids = _text_values(var[idc], f".var[{idc!r}]")
    symc = _pick(var.columns, SYMBOL_COLUMNS, symbol_column)
    symbols = _text_values(var[symc], f".var[{symc!r}]") if symc is not None else keys
    return Axis(list(zip(ids, symbols, keys)), symbols_known=symc is not None)


def read_h5ad_axis(path, id_column: str | None = None, symbol_column: str | None = None,
                   mixed_modalities: str = "refuse") -> Axis:
    try:
        import anndata
    except ImportError as exc:  # pragma: no cover - depends on environment
        raise ValueError("Reading .h5ad requires anndata. Install genejoin with its anndata extra.") from exc
    adata = anndata.read_h5ad(path, backed="r")
    try:
        return axis_from_anndata(adata, id_column, symbol_column, mixed_modalities)
    finally:
        if getattr(adata, "file", None) is not None:
            adata.file.close()


def resolve_input(path) -> Path:
    """Resolve a 10x directory to its unambiguous feature table, or return the file path."""
    path = Path(path)
    if not path.is_dir():
        return path
    directories = (path, path / "filtered_feature_bc_matrix", path / "outs" / "filtered_feature_bc_matrix")
    names = ("features.tsv.gz", "features.tsv", "genes.tsv.gz", "genes.tsv")
    candidates = [folder / name for folder in directories for name in names if (folder / name).is_file()]
    if not candidates:
        raise ValueError(f"{path}: no features.tsv(.gz) or genes.tsv(.gz) found")
    if len(candidates) > 1:
        raise ValueError(f"{path}: multiple 10x feature tables found; supply a specific file")
    return candidates[0]


def load_axis(path, fmt: str = "auto", separator: str = "-", mixed_modalities: str = "refuse",
              id_column: str | None = None, symbol_column: str | None = None) -> Axis:
    """Load one side of a join.

    ``fmt``: ``auto``, ``features`` (raw table; keys generated in file order),
    ``keyed`` (actual final keys supplied) or ``h5ad``.
    """
    path = resolve_input(path)
    if fmt == "auto":
        if path.suffix.lower() == ".h5ad":
            fmt = "h5ad"
        else:
            with _open(path) as handle:
                first = handle.readline().rstrip("\r\n").split("\t")
            fmt = "keyed" if first == ["stable_id", "symbol", "final_key"] else "features"
    if fmt == "h5ad":
        return read_h5ad_axis(path, id_column, symbol_column, mixed_modalities)
    if fmt == "keyed":
        return Axis(read_keyed_table(path))
    if fmt == "features":
        rows = read_feature_table(path, mixed_modalities)
        return Axis([(i, s, None) for i, s in rows], separator)
    raise ValueError(f"Unknown format {fmt!r}")


def name_rows(rows: Sequence[tuple[str, str]], separator: str = "-") -> list[tuple[str, str, str]]:
    """Return ``(id, symbol, key)`` with AnnData-style keys, for export."""
    return [(i, s, k) for (i, s), k in zip(rows, assign_keys([s for _, s in rows], separator))]
