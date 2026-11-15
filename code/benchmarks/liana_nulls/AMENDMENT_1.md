# Amendment 1: a third spatial truth and downstream-response truths

Added on 1 November 2026, after `PLAN.md`. This is an **addition** to the pre-registered plan, not
part of it. Result files written by the added parts record the SHA-256 of this file next to the
plan's.

## What had been seen when this was written

- **Seen:** the T0 fidelity results of this pipeline and the T1 calibration summary, which was
  reused from earlier work.
- **Not yet computed by this pipeline:** any T3 result. The biology stage had not run.
- **Seen in earlier, separate work:** the Open Problems cytokine truth on TNBC. It has 812 rows but
  only about 40 that the tools score, which is why a stronger downstream truth is added here.

## Why

`PLAN.md` judges biology on two spatial truths (TNBC Visium, Allen brain) and the small cytokine
truth. That is thin:
- colocalisation is necessary but not sufficient for signalling;
- one of the two spatial truths has 26 positive pairs.

This amendment adds evidence of two different kinds.

## A1. Third spatial truth, at single-cell resolution: mouse hypothalamic preoptic region

**Source.** Moffitt et al., Science 2018 (doi:10.1126/science.aau5324). Same study, same region,
two assays:
- **scRNA-seq**: GEO GSE113576, about 31,000 cells, 10x. The nulls run on this.
- **MERFISH**: animal 1, 12 slices, 73,655 cells, 161 genes, with the study's own `Cell_class`
  (scverse example data, sha256 371723d4…). Used only as truth.

**Cell classes.** The MERFISH classes, with subclasses merged:
- OD Mature 1–4 → OD Mature;
- OD Immature 1–2 → OD Immature;
- Endothelial 1–3 → Endothelial;
- "Ambiguous" is dropped.

That leaves 9 classes: Inhibitory, Excitatory, Astrocyte, OD Mature, OD Immature, Endothelial,
Ependymal, Microglia, Pericytes.

**scRNA labels by transfer.** The study's per-cell scRNA labels are not downloadable here (the
journal and PMC supplements sit behind a bot check; the CNGB mirror would take about 9 h). So the
scRNA cells get MERFISH classes by transfer:
1. log CP10k; 2,000 highly variable genes; 50 PCs; k-means with k = 40 (seed 0).
2. Each cluster's mean over the genes shared with the MERFISH panel, and each MERFISH class's mean,
   are z-scored per gene across clusters or classes.
3. A cluster takes the class with the highest Pearson correlation, if r ≥ 0.3. Otherwise its cells
   are left out.
4. Canonical markers per class are reported as a check: Gad1, Slc17a6, Aqp4, Mbp, Pdgfra, Cldn5,
   Ccdc153, Cx3cr1, Kcnj8.

The nulls run on the transferred classes with ≥ 20 cells.

**Truth.** Neighbourhood enrichment, per slice:
- graph: 10 nearest neighbours by centroid;
- statistic: number of edges between class s and class t;
- z-score against 200 permutations of the labels within the slice;
- the per-slice z values are averaged over the slices where both classes have ≥ 10 cells;
- colocalised if the mean z ≥ 1.645, the same threshold as the TNBC truth.

**Endpoint.** The primary endpoint of `PLAN.md`, unchanged: the kept-vs-dropped colocalisation odds
ratio among Null A's calls with s ≠ t, with a bootstrap over (s, t) pairs. Plus the secondary
(s, t)-level AUPRC.

## A2. Downstream-response truths (see L → t through the receiver's transcription), all datasets

A call (s, t, L, R) says that t receives L. If t does, t should show the transcriptional footprint
of L's signalling, measured on genes other than L and R. Two footprint resources:
- **CytoSig** (Jiang et al., Nat Methods 2021): `signature.centroid`, 43 cytokine response
  signatures, human. A row is in scope if its ligand complex contains the cytokine's gene. The
  mapping is fixed in `downstream.CYTOKINE_GENES`; NO (nitric oxide) is not a ligand gene and is
  skipped.
- **PROGENy** (Schubert et al., Nat Commun 2018), from OmniPath: 14 pathways, the top 500 genes per
  pathway by p-value (decoupler's convention). An L–R pair is assigned to a pathway, with a sign,
  when every L and R subunit has a PROGENy weight of the same sign there. This is LIANA+'s
  `generate_lr_geneset` rule.
- Mouse datasets use HCOP one-to-one orthologs.

**Activity of a footprint in a cell type t.**
1. The contrast profile of t: mean log-expression in t minus the mean over the other types, per gene.
2. Regress that profile on the footprint weights (univariate linear model, decoupler's `ulm`); the
   activity is the slope's t-value.
3. The genes of L and R are removed from the footprint first.
4. For PROGENy, the activity is multiplied by the L–R's sign.

**Truth.** Per footprint (cytokine, or L–R × pathway), the activities are z-scored across the
receiving types. A unit (ligand or L–R, t) responds if z ≥ 1.645. Sensitivity: z ≥ 1.0.

**Endpoints, per dataset × tool × resource.**
- **D1 (main):** the kept-vs-dropped odds ratio of downstream response among Null A's calls with
  s ≠ t, with a bootstrap over truth units (2,000 resamples). Same form as the primary endpoint.
- **D2 (prediction, fixed now):** a downstream truth sees the receiver, not the sender. So:
  - calls Null B drops because the **receptor** side fails (limiting partner = receptor) should
    respond less often than calls it keeps (odds ratio, kept vs receptor-limited dropped, above 1);
  - calls dropped because the **ligand** side fails need not differ.

  A D2 odds ratio CI entirely below 1 counts against Null B.

## How the additions enter the verdict

The pre-registered decision of `PLAN.md` is computed and reported exactly as written, and is not
changed by this amendment. Next to it the report gives an **extended verdict**:
- criterion (2) of `PLAN.md` evaluated over all three spatial truths (TNBC, Allen brain,
  hypothalamus);
- a downstream line:
  - **supports** if some D1 odds-ratio CI lies above 1 and none lies below 1;
  - **against** if any D1 or D2 CI lies entirely below 1;
  - **neutral** otherwise.

The extended verdict is labelled "Amendment 1" wherever it appears.

## Status

The report keeps a status table, one row per component: whether it is pre-registered or Amendment 1,
and whether it is done, running, pending or failed.
