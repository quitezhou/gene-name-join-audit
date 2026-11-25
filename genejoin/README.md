# genejoin

Check the feature identities on both sides of a gene-name join before combining transcriptomic datasets. Symbols can be unique inside each object while identifying different annotation records across objects.

- When two records share a symbol, filtering each sample before assigning unique names can give the bare symbol to different records. The first check counts shared final names that point to different stable IDs.
- When an annotation release changes a symbol, a shared record can acquire two join keys. An intersection can lose the feature and an outer join can split it across columns. The second check counts shared, individually identifiable IDs that carry different final names.

## Install version 0.1.1

The core uses the Python 3.9+ standard library. The optional `anndata` extra enables reading `.h5ad` files.

```sh
python -m pip install "genejoin @ git+https://github.com/quitezhou/gene-name-join-audit.git@v0.1.1#subdirectory=genejoin"
python -m pip install "genejoin[anndata] @ git+https://github.com/quitezhou/gene-name-join-audit.git@v0.1.1#subdirectory=genejoin"
```

From the repository root, run `python -m pip install ./genejoin` or `python -m pip install "./genejoin[anndata]"`. For offline installation of the core from the release wheel, run `python -m pip install --no-index --no-deps genejoin-0.1.1-py3-none-any.whl`.

## Check the join

```sh
genejoin check reference.h5ad spatial.h5ad --events discordances.tsv
genejoin check sampleA/features.tsv.gz sampleB/outs/filtered_feature_bc_matrix
genejoin check left_keys.tsv right_keys.tsv --left-source reference.tsv.gz --json
```

Use the keys the workflow actually joined. A keyed TSV has a header with the three tab-separated columns `stable_id`, `symbol`, and `final_key`. A raw 10x TSV has two or three columns without a header (ID, symbol, optional feature type). For raw tables, genejoin generates AnnData-style keys on the supplied rows in file order. If a workflow names the complete table before filtering, export the retained rows with their inherited keys. A directory input must identify one unambiguous `features.tsv` or `genes.tsv` file, optionally gzip-compressed.

For AnnData files, the actual `var_names` are the join keys. genejoin opens the file in backed read-only mode. Supply `--id-column` or `--symbol-column` for custom `.var` annotations. Auto-detection recognizes common columns such as `gene_ids` and `gene_symbols`. Verify that the selected ID column identifies comparable records in both datasets. When the ID column is absent, genejoin accepts an index consisting entirely of Ensembl gene IDs. If original symbols are missing, the identity checks still run. Supply complete feature tables via `--left-source` and `--right-source` to identify duplicate families and distinguish them from singleton symbol changes.

Check the number of shared names and shared IDs in the summary. Exit codes are `0` for no detected discordance, `1` for at least one discordance, and `2` for an input or output error. A pair with zero overlap has zero comparable features. Shared IDs present multiple times on either side are counted separately and excluded from the ID-to-name check. The events TSV records individual discrepancies and any original-symbol changes that preserved both the key and the ID. Stable IDs are compared literally, including version suffixes. A different ID establishes a record mismatch. Locus identity and downstream effects require separate evidence.

For CELLxGENE AnnData files the gene/spike-in annotation is stored in `feature_biotype`, and `feature_type` describes gene classes such as `protein_coding` and `lncRNA`. genejoin uses the gene/spike-in annotation to identify gene rows. For 10x AnnData objects it uses the assay modality (`Gene Expression`) in `feature_types`, or in `feature_type` for files using the singular 10x column name. When CELLxGENE biotypes and 10x `feature_types` are both present, gene rows must meet both criteria. [CELLxGENE schema](https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/7.0.0/schema.md#feature_biotype)

Mixed gene and non-gene tables are refused by default when feature types are available. For raw tables, use `--mixed-modalities genes-only` if your workflow selects gene expression rows before naming. For AnnData files, the same option keeps existing keys for the selected gene rows. If naming included multiple modalities, a keyed TSV of the actual retained gene rows also preserves the keys.

## Python

```python
import genejoin

report = genejoin.check_anndata(adata_ref, adata_spatial, id_column="gene_ids")
print(report.same_name_different_id, report.same_id_different_name)
print(report.summary)
for event in report.events:
    print(event)

report = genejoin.check("a/features.tsv.gz", "b/features.tsv.gz")
```

## Recover candidate IDs when a processed object retained only names

```sh
genejoin recover original_features.tsv.gz final_names.txt candidates.tsv
genejoin recover original_features.tsv.gz final_names.txt ordered_candidates.tsv --order-preserved
```

Provide the complete source feature table from the matching reference release and one final key per line in its actual order. Recovery assumes literal symbols with standard AnnData hyphen-number suffixes. The TSV contains one row per candidate ID and the candidate count for that key. A single candidate identifies the record under the supplied assumptions. Multiple candidates preserve uncertainty. Use `--order-preserved` only when provenance confirms that the workflow retained the original ordering. The program rejects orderings that contradict the supplied source.

The manuscript archives include [genejoin v0.1.0](https://github.com/quitezhou/gene-name-join-audit/releases/tag/v0.1.0) for exact reproduction of the study analyses.

## Tests and citation

From this package directory, install the test dependencies with `python -m pip install ".[test,anndata]"` and run `python -m pytest`. The tests include synthetic join and recovery examples. Study-data checks use external feature tables where available. See `tests/data/README.md`.

Cite the software version using `CITATION.cff`. The MIT license in this directory covers genejoin.
