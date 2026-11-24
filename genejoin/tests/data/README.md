Optional real-data tests reproduce the manuscript's probe-panel example.

Either copy the two 10x feature tables here as

* `reference_36601.tsv.gz` (36,601-row GRCh38 2020-A reference)
* `ffpe_panel_17943.tsv.gz` (17,943-row FFPE probe panel, an exact ordered subset)

or point the tests at the manuscript archive's `systematic_mechanism/panel_axes` folder:

```shell
GENEJOIN_PANEL_AXES=/path/to/systematic_mechanism/panel_axes pytest
```

The tests select the two named source axes listed in the manuscript's `panel_sources.tsv` and check that ARMCX5-GPRASP2 is the single
same-name, different-ID pair, and that order-constrained recovery leaves only MATR3 ambiguous.

To also run the lung and blood regressions, set `GENEJOIN_STUDY_DIR` to the extracted
`audit_20261123/gb_followup` directory. For the H5AD comparisons, set
`GENEJOIN_PBMC3K_H5AD` and `GENEJOIN_SYMBOL_ONLY_H5AD` to the two downloaded files
specified in `tests/test_study_inputs.py`. The tests check their SHA-256 hashes. The
test suite does not download any data.
