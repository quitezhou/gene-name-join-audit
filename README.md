# Auditing gene-name joins in transcriptomic datasets

Analysis code for *Joining transcriptomic datasets by gene name can silently pair different genes, even from one reference*. The audit checks whether a shared final name points to different source IDs and whether one shared ID appears under different final names.

## Check your own join with genejoin

For every shared name, genejoin checks whether both datasets assign it the same stable ID. For every shared, individually identifiable ID, it checks whether both datasets give it the same name.

Install the [version 0.1.1 release wheel](https://github.com/quitezhou/gene-name-join-audit/releases/tag/v0.1.1) with Python 3.9 or later. The core has no runtime dependencies.

```sh
python -m pip install "https://github.com/quitezhou/gene-name-join-audit/releases/download/v0.1.1/genejoin-0.1.1-py3-none-any.whl"
genejoin check sampleA/features.tsv.gz sampleB/features.tsv.gz
```

For AnnData files install the optional dependency.

```sh
python -m pip install "genejoin[anndata] @ https://github.com/quitezhou/gene-name-join-audit/releases/download/v0.1.1/genejoin-0.1.1-py3-none-any.whl"
genejoin check reference.h5ad spatial.h5ad --events discordances.tsv
```

The checker exits with status 0 for no detected discordance, 1 for a discordance and 2 for an input or output error. Check the overlap counts to see how many records could be compared. The Python interface accepts in-memory objects via `genejoin.check_anndata(adata_a, adata_b)`.

The [package guide](genejoin/README.md) describes recovery when a processed object lacks stable IDs. The [test input guide](genejoin/tests/data/README.md) explains how to reproduce the probe-panel, lung, blood and AnnData checks. The [MIT license](genejoin/LICENSE) covers the genejoin package.

## Try a pair of gene tables

From the repository root, run the checker with Python 3.10 or later. This example needs only the standard library.

```sh
python audit_20261123/gb_followup/systematic_mechanism/audit_join.py examples/full.tsv examples/subset.tsv --events events.tsv
```

In the synthetic full table, both ID_A and ID_B have the symbol GENE_X. AnnData-style local naming assigns GENE_X to ID_A and GENE_X-1 to ID_B. The ordered subset retains only ID_B and assigns it GENE_X. Both names in the smaller table occur in the larger one. The output records one shared name with different IDs and one shared ID with different names. The program prints summary counts and writes both events to `events.tsv`.

For a new analysis, provide two headerless, tab-separated gene tables in their actual row order, with source identifier and raw symbol in the first two columns. A third feature-type column may contain `Gene Expression`. If final keys already exist, use the keyed-input options described by `--help`; keyed files have the header `stable_id<TAB>symbol<TAB>final_key`. Original prefilter tables can establish which symbols were duplicated upstream. Literal ID disagreement establishes a record discrepancy. Distinguishing disjoint loci, historical annotation revisions and effects on a biological endpoint requires further evidence.

## Reproduce the study figures

The paper and its supplementary archives use the frozen [genejoin v0.1.0 release](https://github.com/quitezhou/gene-name-join-audit/releases/tag/v0.1.0). For new analyses use v0.1.1, which also reads CELLxGENE gene and spike-in annotations. The standalone example above works from this repository alone. To reproduce the study analyses, obtain both supplementary archives from the authors or the article when available. Extract **both into the same separate working directory**. This retains the inputs, manifests, archived third-party materials and exact code snapshot used in the analysis. From that directory run:

```sh
python -m pip install -r requirements-review.txt
python offline_check.py
```

The package check uses Python 3.11. It verifies the payload against its SHA-256 manifest, replays the supported analyses and renders Figures 1 and 2 and Supplementary Figure S1. The PNG comparison uses the rendering environment described in the supplements. Fresh fits from original study inputs may also need R, limma, AnnData or PyDESeq2 as specified in the study methods. The archived README describes these optional source-data replays. The checksums below identify the two inputs for this code snapshot.

| Supplement | SHA-256 |
| --- | --- |
| Additional_file_2.zip | `a35c7d700032f6bdb8465242763405636dea67f6193238ee5c46746aff3e5dd0` |
| Additional_file_3.zip | `1603f7dfea361d043632b3159267ac2b4c7e119dcb9f3489a16dc171561e65e0` |

## Code map

| Analysis | Entry point |
| --- | --- |
| Name-to-ID and ID-to-name checks | [audit_join.py](audit_20261123/gb_followup/systematic_mechanism/audit_join.py) |
| Candidate IDs for a table with missing identifiers | [recover_ids.py](audit_20261123/gb_followup/prejoin_bounds/recover_ids.py) |
| Skin and cortex processing order | [skin replay](audit_tangram_skin_real_counts.py), [cortex replay](audit_20261123/gb_followup/external_order_replication/replay.py) |
| Filtering thresholds in eight blood donors | [PBMC replay](audit_20261123/gb_followup/prospective_holdout/verify.py) |
| Coverage and cost of matching strategies | [comparison](audit_20261123/gb_followup/repair_benchmark/benchmark.py) |
| Published adipose marker, prostate treatment contrast and lung colocalization | [adipose bounds](audit_20261123/gb_followup/marker_resolution/bound_contrast.py), [prostate refit](audit_20261123/gb_followup/full_S11/refit.py), [lung replay](audit_20261123/gb_followup/fifth_case/replay_figure3d.py) |
| Manuscript figures | [figure scripts](audit_20261123/gb_followup/v2/submission/figures_final.py) |

Historical dataset versions, source-file digests, accession lists and study-specific data references are documented in the supplements. Public deposits include the [blood donor experiment](https://doi.org/10.1038/nbt.4042), [cortex study](https://doi.org/10.1038/s41593-020-00787-0), [adipose atlas](https://doi.org/10.1038/s41467-023-36983-2), [lung cancer study](https://doi.org/10.1038/s41467-024-48700-8) and [prostate cancer study](https://doi.org/10.1038/s41467-024-54364-1). Original software includes [AnnData](https://github.com/scverse/anndata) and [Tangram](https://github.com/broadinstitute/Tangram). The examples in this repository are synthetic.
