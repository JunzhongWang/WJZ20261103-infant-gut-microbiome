# Convergence without loss of individuality in the developing infant gut microbiome

Analysis code, data tables and result files accompanying the manuscript submitted to *mSystems*.

**Version:** v1.1.0 · **Released:** 2026-11-13 · **Contents:** 189 files, 24 MB
**Code:** Python 3 (standard library only) + two base-R helpers · **License:** MIT (code), CC BY 4.0 (derived tables and data matrices)

This archive lets a reader go from public data to every number, figure and table in the paper
without contacting the authors, and it includes the analysis-ready data themselves: the seven
per-cohort species-level abundance matrices and their sample metadata, 3,142 samples from 947
infants, in `data/<cohort>/` (section 2.2). Raw sequence reads are **not** redistributed —
they stay in the public archives under the BioProject accessions listed below.

---

## 1. Repository layout

| Path | What it is |
|---|---|
| `pipeline/` | Final analysis pipeline. Flat directory on purpose: several scripts import each other as modules (`decay_model`, `make_figs_v3`, `fetch_cmd3_direct`, `infant_identification`), so splitting them into sub-packages breaks the import chain |
| `pipeline/archive/` | Superseded iterations and the exploratory MGnify screening path, kept for provenance only. Not needed to reproduce the paper |
| `results/real/` | Every numerical result of the core age/interval analysis (JSON/CSV), all run logs, and the figure files produced by the plotting scripts |
| `results/real/png/` | Raster versions of the same panels |
| `results/ibd2/` | The two adult analysis tables (`cohort_*.csv` + `abund_*.csv` for HMP2 and Hall 2017) that Figures 6–7 need, prepared by the companion adult analysis and shipped verbatim so the infant–adult contrast comes from one estimator |
| `results/comparison/` | The infant-vs-adult R contrast (Figure 6): per-cohort and pooled R by interval band, the infant/adult ratio with CI, and the three sensitivity checks |
| `results/identification/` | The re-identification experiment (Figure 7): recognition rates per band and group, full-pool sensitivity, and the complete trial-level dumps |
| `results/decomposition/` | The taxonomic decomposition of the within-infant distance (Figure 8): three components per band, per-genus contributions and deltas, acquisition/loss event rates |
| `tables/` | Source CSV for Table 1, Table 2, Table 3 and the taxon-level supplementary table |
| `data/<cohort>/` | **Analysis-ready data:** per-cohort `community.csv` (leaf partition, samples × species relative abundance), `species.csv` (all MetaPhlAn ranks, used by the species-rank sensitivity checks) and `meta.csv` (age, subject, antibiotic flag, country), at the exact paths `pipeline/` reads. `data/INDEX.csv` lists dimensions, accessions and SHA-256; `data/README.md` documents the columns |
| `figures_tables_standalone/` | Manuscript-ready exports: Figures 1–8 (SVG + PNG), Figure S1 (SVG/PNG/PDF + caption), Tables 1–3 (HTML + CSV) |
| `docs/exploratory/` | Early cohort-screening report (Chinese), retained as an audit trail |
| `MANIFEST.csv` | Every file with byte size and SHA-256, for integrity checking |
| `CITATION.cff` | Machine-readable citation metadata |
| `CHANGELOG.md` | What changed between deposit versions |

No file contains an absolute path from the authoring machine: every script resolves `data/`,
`results/` and `figures_tables_standalone/` relative to the directory that contains `pipeline/`.
Run everything from the archive root and nothing needs editing.

---

## 2. Data sources

### 2.1 Primary data: curatedMetagenomicData v3

All taxonomic and functional profiles come from **curatedMetagenomicData v3** (cMD3,
Bioconductor ExperimentHub), which uniformly reprocesses public infant cohorts with the
**bioBakery 3** pipeline (MetaPhlAn 3 for taxonomy, HUMAnN 3 for pathways). One profiler for
every cohort removes the largest source of cross-study technical variation.

Seven infant stool cohorts are used. ExperimentHub resource IDs, release dates and the original
BioProject accessions are:

| Cohort | BioProject (raw reads) | ExperimentHub ID | cMD3 release | Country | Design |
|---|---|---|---|---|---|
| Shao 2019 | **PRJEB32631** | EH5800 | 2021-03-31 | GBR | sparse, term infants |
| Vatanen 2016 (DIABIMMUNE) | **PRJNA290380** | EH5860 | 2021-03-31 | RUS | sparse, prospective birth cohort |
| Bäckhed 2015 | **PRJEB6456** | EH7097 | 2021-10-14 | SWE | sparse, mother–infant |
| Yassour 2018 | **PRJNA475246** | EH5914 | 2021-03-31 | FIN | sparse, family cohort |
| Kostic 2015 (TEDDY) | **PRJNA231909** | EH5638 | 2021-03-31 | FIN/USA | dense pilot |
| Ferretti 2018 | **PRJNA352475** | EH5536 | 2021-03-31 | ITA | dense, every 4 days |
| Brooks 2017 | **PRJNA376566** | EH7109 | 2021-10-14 | USA | dense preterm NICU |

Functional profiles (Supplementary Figure S1) are the matching `<study>.pathway_abundance`
resources from the same releases, retrieved by `pipeline/fetch_pathways.py`.

**Retrieval route.** `pipeline/fetch_cmd3_direct.py` speaks to ExperimentHub directly: it
downloads the hub catalogue (`experimenthub.sqlite3`, ~10 MB), selects by SQL the newest resource
whose title matches `<study>.relative_abundance`, downloads that `.rda` with HTTP-range resume,
and hands it to `pipeline/rda_extract.R`. No Bioconductor client, no R package installation. Two
releases are in play because the script always takes the newest available resource: five cohorts
resolve to the 2021-03-31 build, Bäckhed 2015 and Brooks 2017 to 2021-10-14.

**Processing already done upstream.** Profiling and quality control are those of the cMD3
pipeline; raw reads were not reprocessed here. Everything in `pipeline/` operates on the
processed abundance tables.

**Resulting analysis set.** 3,842 retrieved metagenomes → 96 non-stool samples, 41 samples
outside the 0–1,095-day window and 563 without a usable age removed → **3,142 stool samples from
947 infants**, ages 4–1,093 days. Ages for 356 samples were reconstructed from the subject's
median offset between `infant_age` and `days_from_first_collection` (353 of them in Vatanen 2016,
which reports `days_from_first_collection` for all 785 of its stool samples but `infant_age` for
only 293).

### 2.2 Data included in this archive

`data/<cohort>/` holds the analysis-ready tables, one directory per cohort, and they are the
exact files every script in `pipeline/` reads — byte-identical to the working copies, not
re-derived for the deposit:

| Cohort | `community.csv` | `meta.csv` |
|---|---|---|
| BackhedF_2015 | 287 samples × 533 species | 287 rows |
| BrooksB_2017 | 408 × 225 | 408 |
| FerrettiP_2018 | 97 × 436 | 97 |
| KosticAD_2015 | 118 × 332 | 118 |
| ShaoY_2019 | 1,470 × 687 | 1,470 |
| VatanenT_2016 | 607 × 585 | 607 |
| YassourM_2018 | 155 × 244 | 155 |

Together 12 MB, 3,142 samples, 947 infants. Each cohort directory now also carries
`species.csv`, the full MetaPhlAn rank hierarchy (kingdom through species) that the
species-rank sensitivity checks of Sections 2.7 and 2.8 read; `data/INDEX.csv` records the
per-file dimensions, the BioProject accession, the ExperimentHub resource ID and the SHA-256
of each file; `data/README.md` documents the column semantics. These tables are released
under CC BY 4.0.

### 2.3 Adult cohorts used for the contrast (Figures 6–7)

Two adult cohorts provide the reference side of the infant–adult comparison: the IBD
Multi'omics deposit of the Human Microbiome Project 2 (BioProject **PRJNA398089**) and the
cohort of Hall 2017 (**PRJNA385949**), retrieved through the same ExperimentHub resources as
the infant tables. Their analysis-ready tables ship in `results/ibd2/`
(`cohort_<study>.csv` with sample metadata and group labels, `abund_<study>.csv` with the
species-rank abundance matrix) at exactly the paths `infant_vs_adult_r.py` and
`infant_identification.py` read. They were prepared with the companion adult analysis
(deposited separately) and are included verbatim here so that both sides of every
infant–adult contrast come from one estimator, one distance and one pipeline. Only the
non-IBD (HMP2) and healthy (Hall) subsets enter the figures; the CD/UC/IBD subsets are in
the tables for completeness.

What is *not* included is the 400 MB download cache (`.rda` blobs from ExperimentHub plus
intermediate TSVs); `pipeline/fetch_cmd3_direct.py` recreates it, and Stage 1 of §5 does so
end-to-end. The retrieval is only needed if a reader wants to rebuild the tables from source
rather than start from `data/<cohort>/`.

---

## 3. Computational environment

| Item | Requirement | Tested on |
|---|---|---|
| Python | **≥ 3.9**, standard library only | CPython 3.13.12 (macOS arm64) |
| Third-party Python packages | **none** — no numpy, pandas, scipy, matplotlib, biopython | — |
| R | **≥ 4.0**, base R only (no Bioconductor packages) | R 4.1.0 |
| Disk | ~400 MB for the download cache | — |
| Network | HTTPS to `experimenthub.bioconductor.org` (Azure blob storage behind it) | — |
| Hardware | any laptop; no GPU | — |

Standard-library modules used across the pipeline: `argparse`, `csv`, `glob`, `importlib`,
`json`, `math`, `os`, `pathlib`, `random`, `re`, `shutil`, `sqlite3`, `statistics`,
`subprocess`, `sys`, `time`, `urllib`.

Everything quantitative is implemented from scratch on the standard library: the Bray–Curtis,
Jaccard and Aitchison distance matrices, the sparse-matrix shortcuts used for the pathway tables,
the bootstrap resampling, the DerSimonian–Laird random-effects model, and the SVG renderer. That
is deliberate — a reviewer can read the whole calculation without installing a scientific stack.

R is needed only for the two `.rda` extraction helpers (`pipeline/rda_extract.R`,
`pipeline/rda_extract_pathway.R`), which use `load()` and `attr()` on the plain matrix or
`ExpressionSet` that cMD3 stores. If you start from the shipped `results/`, R is not needed at
all.

Runtime: the download dominates (~10–20 min on a domestic connection, resumable and idempotent);
all computation afterwards is a few minutes. Bootstrap counts are command-line arguments: 1,000
for the decay and strata models, 2,000 for the meta-analysis, 500 for the functional layer.

---

## 4. Script reference

Every script is run from the archive root. `ROOT` inside each script is the directory holding
`pipeline/`, so relative paths resolve without configuration.

### 4.1 Retrieval (stage 1)

| Script | Purpose | Inputs | Outputs |
|---|---|---|---|
| `fetch_cmd3_direct.py` | **Primary retrieval.** ExperimentHub catalogue → newest `<study>.relative_abundance` → resumable `.rda` download → calls `rda_extract.R`. Idempotent; caches under `data/cache/` | ExperimentHub (HTTPS) | `data/cache/*.rda`, `data/cache/sample_metadata.csv`, `data/<study>/{species.tsv, pheno.tsv}` |
| `rda_extract.R` | Turn a `relative_abundance` `.rda` into a mass-preserving **leaf-only** table (percent, rows sum to 100) | `<study>` `.rda` | `data/<study>/{species.tsv, pheno.tsv}` |
| `fetch_pathways.py` | Same retrieval for the `<study>.pathway_abundance` resources used by the functional layer | ExperimentHub | `data/<study>/pathway.tsv` |
| `rda_extract_pathway.R` | Turn a `pathway_abundance` `.rda` into a sample × pathway TSV | `<study>` pathway `.rda` | `data/<study>/pathway.tsv` |
| `check_network.sh` | Diagnostic: which ExperimentHub hosts this machine can reach, and what to do if not | — | console report |
| `fetch_cmd3.R` | *Alternative*, not used for the submitted results: retrieval through the `curatedMetagenomicData` R package (requires a Bioconductor installation) | ExperimentHub via R | `data/<study>/species.tsv` |

### 4.2 Cohort construction (stage 2)

| Script | Purpose | Inputs | Outputs |
|---|---|---|---|
| `build_cohort.py` | Build the analysis tables: body-site filter, 0–1,095-day age window, age re-anchoring within subject, percent → proportion | `data/<study>/species.tsv`, `data/cache/sample_metadata.csv` | `data/<study>/{community.csv, meta.csv}` |
| `rebuild_meta.py` | Rebuild `meta.csv` on the birth-anchored `infant_age` axis; keeps `days_from_first_collection` alongside as the sensitivity axis | `data/cache/sample_metadata.csv` | `data/<study>/meta.csv` |
| `join_meta.py` | *Exploratory path*: join an MGnify species table with cMD3 metadata by run accession | MGnify `species.tsv`, `infant_samples.csv` | `data/<study>/{species.csv, meta.csv}` |

### 4.3 Statistical models (stage 3)

| Script | Purpose | Inputs | Outputs |
|---|---|---|---|
| `decay_model.py` | Pair-level regression separating **age** from **sampling interval**: within- and between-subject Bray–Curtis distance as a function of both, with bootstrap CIs (Figure 1) | `community.csv`, `meta.csv` | `decay_model.{csv,json,log}` |
| `strata.py` | **Model-free** estimate of R (between/within distance ratio) in cells matched on age *and* interval; cell means with bootstrapped intervals (Table 2) | `community.csv`, `meta.csv` | `strata.{csv,json,log}` |
| `meta_decay.py` | DerSimonian–Laird random-effects meta-analysis of per-cohort age slopes; heterogeneity I² (Figure 2) | per-cohort statistics | `meta_slopes.csv`, `meta_buckets.csv`, `meta_decay.{json,log}` |
| `joint_model.py` | Joint model of age and interval over every sample pair; within- vs between-subject response surface | `community.csv`, `meta.csv` | `joint_curve.csv`, `joint_surface.csv`, `joint_model.{json,log}` |
| `identify.py` | Design-based identification check needing no extrapolation: dense cohorts answer the age/interval question by design alone | `community.csv`, `meta.csv` | `identify.{csv,json,log}` |
| `taxon_icc.py` | Per-taxon individual signature: variance partition (ICC) on CLR abundances, one cohort at a time (Figure 5, Table S1) | `community.csv` | `taxon_icc.tsv` (per cohort, redirected to `<study>.taxon.txt`) |
| `revision.py` | Metric sensitivity: repeat the whole comparison with Bray–Curtis, Jaccard (presence–absence) and Aitchison (CLR-Euclidean) (Figure 4) | `community.csv`, `meta.csv` | `revision.{json,log}` |
| `revision_trends.py` | Age trends of within-subject distance, between-subject distance and R separately, inside matched strata, inverse-variance weighted (Figure 3, Table 3) | `strata.json` | `age_trends_primary.{csv,json}`, `age_trends_term_only.{csv,json}` |
| `functional_R.py` | Functional-layer robustness: recompute R on HUMAnN pathway abundances, 500 bootstraps over infants (Figure S1) | `pathway.tsv`, `meta.csv` | `functional_R.json` |
| `infant_pipeline.py` | Original single-cohort summary (`analyze --dir`). Superseded by the model scripts above, but it produces the per-cohort JSON shipped in `results/real/` | `community.csv`, `meta.csv` | `results/real/<study>.json` |

### 4.4 Figures and tables (stage 4)

| Script | Purpose | Inputs | Outputs |
|---|---|---|---|
| `make_figs_v3.py` | SVG panel library: sampling design, R against interval, species-level ICC (680 px vector panels, English labels) | model outputs | `fig1_cohorts.svg`, `fig2_ratio.svg`, `fig4_icc.svg` |
| `make_figs_v4.py` | Revision figures: the two distances falling together with age, interval dependence of R inside each cohort, ICC summary | `strata.json`, `revision.json`, `age_trends_primary.json` | `fig3_agetrends.svg`, `fig4_cohortmetric.svg`, `fig5_icc.svg` |
| `make_figS1.py` | Render Supplementary Figure S1 (functional robustness) | `functional_R.json` | `figS1_functional.svg` |
| `tables.py` | Table 1 (cohort characteristics) and the species × cohort ICC matrix | `community.csv`, `meta.csv`, `taxon_icc.tsv` | `table1.csv`, `table_taxon.csv` |
| `export_standalone.py` | Copy figures and rebuild the three tables as manuscript-ready files | `results/real/*` | `figures_tables_standalone/` (Figures 1–5, Tables 1–3) |
| `export_figS1.py` | Export Figure S1 as a deliverable set | `figS1_functional.svg`, `functional_R.json` | `FigureS1_functional_robustness.{svg,png,pdf,html}`, `_caption.txt` |
| `run_all.sh` | Original entry point: fetch + per-cohort summaries for four cohorts. A smoke test, not the full run | — | `results/real/<study>.json` |

### 4.5 Archived scripts

`pipeline/archive/` holds `make_figures.py`, `make_fig_v2.py`, `make_fig_real.py` (earlier figure
generations, superseded by `make_figs_v4.py`) and `mgnify_find.py`, `mgnify_screen.py`,
`mgnify_fetch.py` (the initial MGnify route to candidate cohorts, abandoned in favour of cMD3).
They are kept so the cohort choice is auditable; they are not part of the reproduction path, and
their `ROOT` resolution assumes they sit one level higher.

### 4.6 The three additions of v1.1.0 (adult contrast, re-identification, decomposition)

These six scripts extend the core analysis with the comparisons that date the emergence and
name the carriers of individuality. All six read the same `data/<cohort>/` tables (and, for
the adult side, `results/ibd2/`) and write into their own `results/` subdirectory.

| Script | Purpose | Inputs | Outputs |
|---|---|---|---|
| `infant_vs_adult_r.py` | **Adult contrast (Figure 6).** The matched-interval R estimator applied unchanged to the seven infant cohorts and recomputed for the two adult cohorts; pooled infant vs adult R per band, the infant/adult ratio with CI, and three sensitivity checks (older-sample age stratum, species-rank partition, term-only) | `data/<cohort>/`, `results/ibd2/` | `results/comparison/` (7 files; `contrast_infant_adult.csv` is the key output) |
| `fig_infant_vs_adult.py` | Render Figure 6 (two panels: R against interval, infants vs adults) | `results/comparison/` | `fig_infant_vs_adult.svg`/`.png` |
| `infant_identification.py` | **Re-identification experiment (Figure 7).** 1-nearest-neighbour match of an earlier sample against ten age-matched candidates (chance exactly 10%), per interval band, infants vs adults; fixed-pool primary read-out and full-natural-pool sensitivity; complete per-trial dumps | `data/<cohort>/`, `results/ibd2/` | `results/identification/` (9 files; `ident_groups.csv` is the key output) |
| `fig_identification.py` | Render Figure 7 (two panels: per-cohort curves, infants vs adults) | `results/identification/` | `fig_identification.svg`/`.png` |
| `infant_taxonomic_decomposition.py` | **Taxonomic decomposition (Figure 8).** Exact per-taxon split of the within-infant Bray–Curtis distance into colonisation / extinction / shift of residents, per genus and interval band, with acquisition/loss event rates and the age-stratified sensitivity | `data/<cohort>/` | `results/decomposition/` (10 files; `components_by_band.csv`, `taxa_delta.csv` are the key outputs) |
| `fig_taxonomic_decomposition.py` | Render Figure 8 (three panels: stacked components, contribution deltas, event rates) | `results/decomposition/` | `fig_taxonomic_decomposition.svg`/`.png` |

Monte-Carlo note: every estimator in this archive resamples (cluster bootstrap, 15,000
sampled between-infant pairs, and the 9-other draw of the fixed candidate pool). The shipped
result files are the realisation reported in the manuscript; re-running the pipeline
reproduces every estimate to within about 1% relative on R and 3 percentage points on
recognition rates, with confidence intervals of the same width. The identification summary
tables can be recomputed exactly from the shipped trial dumps
(`ident_trials_fixed10.csv.gz`), which record every trial's pool, rank and hit. Since v1.1.0
the fixed-pool draw is seeded from the trial key (`crc32`), so it no longer depends on the
process hash seed and is identical across machines and runs.

---

## 5. Run order

Stages 1–2 must run in that order. Within stages 3–4 the order matters only where one script's
output is another's input (`strata.py` → `revision_trends.py`; `functional_R.py` →
`make_figS1.py`).

```bash
# ---- stage 1: retrieval (~10-20 min, ~400 MB; resumable, idempotent) -------------
python3 pipeline/fetch_cmd3_direct.py --study ShaoY_2019 VatanenT_2016 BackhedF_2015 \
        YassourM_2018 KosticAD_2015 FerrettiP_2018 BrooksB_2017
#   -> data/cache/{*.rda, sample_metadata.csv}, data/<study>/{species.tsv, pheno.tsv}
python3 pipeline/fetch_pathways.py          # needed only for Figure S1

# ---- stage 2: cohort construction (< 1 min) ------------------------------------
python3 pipeline/build_cohort.py --study ShaoY_2019 VatanenT_2016 BackhedF_2015 \
        YassourM_2018 KosticAD_2015 FerrettiP_2018 BrooksB_2017
#   -> data/<study>/{community.csv, meta.csv}     (3,842 -> 3,142 samples)

# ---- stage 3: statistics (a few minutes total) ---------------------------------
python3 pipeline/identify.py            # -> identify.{csv,json}
python3 pipeline/decay_model.py         # -> decay_model.{csv,json}
python3 pipeline/meta_decay.py          # -> meta_slopes.csv, meta_decay.json
python3 pipeline/strata.py              # -> strata.{csv,json}         [before revision_trends]
python3 pipeline/revision_trends.py     # -> age_trends_primary.{csv,json}
python3 pipeline/joint_model.py         # -> joint_curve.csv, joint_surface.csv
python3 pipeline/revision.py            # -> revision.json
for s in ShaoY_2019 VatanenT_2016 BackhedF_2015 YassourM_2018 \
         KosticAD_2015 FerrettiP_2018 BrooksB_2017; do
  python3 pipeline/taxon_icc.py --dir data/$s > results/real/$s.taxon.txt
done
python3 pipeline/functional_R.py        # -> functional_R.json        [needs stage 1b]

# ---- stage 3b: the three additions of v1.1.0 (~2 min total) --------------------
python3 pipeline/infant_vs_adult_r.py --boot 500   # -> results/comparison/
python3 pipeline/infant_identification.py --boot 1000   # -> results/identification/
python3 pipeline/infant_taxonomic_decomposition.py --boot 500   # -> results/decomposition/

# ---- stage 4: figures, tables, exports (< 1 min) ------------------------------
python3 pipeline/make_figs_v3.py        # -> fig1_*, fig2_ratio.svg, fig4_icc.svg
python3 pipeline/make_figs_v4.py        # -> fig3_agetrends.svg, fig4_cohortmetric.svg, fig5_icc.svg
python3 pipeline/make_figS1.py          # -> figS1_functional.svg
python3 pipeline/fig_infant_vs_adult.py        # -> results/comparison/fig_infant_vs_adult.svg
python3 pipeline/fig_identification.py         # -> results/identification/fig_identification.svg
python3 pipeline/fig_taxonomic_decomposition.py  # -> results/decomposition/fig_taxonomic_decomposition.svg
python3 pipeline/tables.py              # -> table1.csv, table_taxon.csv
python3 pipeline/export_standalone.py   # -> figures_tables_standalone/ (Figures 1-8, Tables 1-3)
python3 pipeline/export_figS1.py        # -> figures_tables_standalone/FigureS1_*
```

Bootstrap count is set per script with `--boot` (defaults: 1,000 for decay and strata, 2,000 for
the meta-analysis, 500 for the functional layer). Reseeding a run moves the third decimal, not
the conclusions.

### Output → manuscript map

| Manuscript item | Producing script | Key result file |
|---|---|---|
| Figure 1 — sampling design and apparent decay | `decay_model.py` → `make_figs_v3.py` | `decay_model.csv`, `fig1_cohorts.svg` |
| Figure 2 — R against sampling interval | `meta_decay.py` → `make_figs_v3.py` | `meta_slopes.csv`, `fig2_ratio.svg` |
| Figure 3 — the two distances fall together with age | `revision_trends.py` → `make_figs_v4.py` | `age_trends_primary.json`, `fig3_agetrends.svg` |
| Figure 4 — interval dependence inside each cohort, three metrics | `revision.py` → `make_figs_v4.py` | `revision.json`, `fig4_cohortmetric.svg` |
| Figure 5 — species-level individual signature | `taxon_icc.py` → `make_figs_v4.py` | `<study>.taxon.txt`, `fig5_icc.svg` |
| Figure 6 — the adult signature is absent in infancy | `infant_vs_adult_r.py` → `fig_infant_vs_adult.py` | `results/comparison/contrast_infant_adult.csv`, `fig_infant_vs_adult.svg` |
| Figure 7 — re-identification accuracy against interval | `infant_identification.py` → `fig_identification.py` | `results/identification/ident_groups.csv`, `fig_identification.svg` |
| Figure 8 — taxonomic decomposition of the within-infant distance | `infant_taxonomic_decomposition.py` → `fig_taxonomic_decomposition.py` | `results/decomposition/components_by_band.csv`, `taxa_delta.csv` |
| Figure S1 — functional (pathway) robustness | `functional_R.py` → `make_figS1.py` | `functional_R.json`, `figS1_functional.svg` |
| Table 1 — cohort characteristics | `identify.py` → `tables.py` | `table1.csv` → `tables/Table1_cohorts.csv` |
| Table 2 — R in matched strata | `strata.py` → `export_standalone.py` | `strata.csv` → `tables/Table2_strata_R.csv` |
| Table 3 — age slopes of the three quantities | `revision_trends.py` → `export_standalone.py` | `age_trends_primary.csv` → `tables/Table3_age_slopes.csv` |
| Table S1 — taxon-level ICC matrix | `taxon_icc.py` → `tables.py` | `table_taxon.csv` → `tables/TableS1_taxon_level.csv` |

---

## 6. Implementation notes

Three choices explain most of the difference between this pipeline and a default workflow. Each
is documented at the top of the script that implements it.

* **Leaf taxa only.** cMD3 stores every rank of the MetaPhlAn tree, so each sample carries seven
  nested copies of the same community and the columns sum to 700, not 100. `rda_extract.R` keeps
  only rows without descendants — those, and only those, form a true partition. Distances
  computed on the full table would double-count ancestors and over-weight coarse ranks.
* **No rarefaction.** Relative abundances are used as provided. Rarefying compositional data
  adds variance without removing the bias it is meant to address.
* **Two age axes.** `infant_age` (days since birth) is primary; `days_from_first_collection` is
  carried as a sensitivity axis (`revision_trends.py --tag term_only`). The original pass used
  the latter by mistake — in Vatanen 2016 the median gap between the two is 139 days, so that run
  was plotting time-since-enrolment against a birth-anchored axis. The comparison is reported
  rather than quietly dropped.

Note that `tables.py`, `make_figs_v3.py`, `make_figs_v4.py`, `make_figS1.py`,
`export_standalone.py` and `export_figS1.py` execute at import time (no `main()` guard) and
therefore require `data/` and `results/real/` to exist. They are scripts, not importable
libraries.

---

## 7. What is not in this archive

* `data/` — regenerated by stage 1 (~400 MB of `.rda` cache plus derived TSVs), and listed in
  `.gitignore`.
* Manuscript typesetting scripts (citation renumbering, HTML→DOCX/PDF export, figure embedding).
  They operate on the submitted HTML source and are not part of the analysis.

---

## 8. License and citation

Code is released under the MIT license (`LICENSE`); derived tables and result files under
CC BY 4.0. If you use this archive, please cite the article and this record — metadata are in
`CITATION.cff`.

**DOI:**&#8203; [10.5281/zenodo.23906301](https://doi.org/10.5281/zenodo.23906301) (concept DOI, resolves to the latest version; this release: [10.5281/zenodo.24187660](https://doi.org/10.5281/zenodo.24187660))

**Code repository:** [github.com/JunzhongWang/WJZ20261103-infant-gut-microbiome](https://github.com/JunzhongWang/WJZ20261103-infant-gut-microbiome)
