"""genejoin: check gene-name joins between transcriptomic datasets before you trust them.

Two checks, run before joining two datasets by gene name:

1. every shared name should denote the same stable ID on both sides, and
2. every shared stable ID should carry the same name on both sides.

Typical use::

    import genejoin
    report = genejoin.check_anndata(adata_a, adata_b)
    if not report.ok:
        print(report)
        for event in report.events: ...
"""

from __future__ import annotations

from .core import Axis, Reference, Report, assign_keys, audit, recover
from .io import axis_from_anndata, load_axis, read_feature_table, read_keyed_table, read_keys

__version__ = "0.1.1"
__all__ = ["Axis", "Reference", "Report", "assign_keys", "audit", "recover", "check", "check_anndata",
           "axis_from_anndata", "load_axis", "read_feature_table", "read_keyed_table", "read_keys"]


def check(left, right, **kwargs) -> Report:
    """Audit two inputs given as paths (feature tables, keyed tables, 10x folders or .h5ad)."""
    return audit(load_axis(left, **kwargs), load_axis(right, **kwargs))


def check_anndata(left, right, id_column: str | None = None, symbol_column: str | None = None,
                  mixed_modalities: str = "refuse") -> Report:
    """Audit two in-memory AnnData objects using ``var_names`` and their ID columns."""
    return audit(axis_from_anndata(left, id_column, symbol_column, mixed_modalities),
                 axis_from_anndata(right, id_column, symbol_column, mixed_modalities))
