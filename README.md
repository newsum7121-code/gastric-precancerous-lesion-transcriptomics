# Evidence-stratified transcriptomic prioritization in gastric precancerous lesions

This repository contains a limited, versioned set of analysis scripts and aggregate gene-level outputs for the manuscript *Evidence-stratified prioritization of transcriptomic candidates across gastric precancerous lesions and cancer: a cross-cohort public-data analysis*.

## Scope

The included scripts reproduce the GSE163416 gene-level edgeR sensitivity refit, the five-candidate leave-one-sample-out analysis, and the GO Biological Process sensitivity analysis, using public processed inputs and the historical comparator table described below. These are the modules locally rerun against the manuscript-linked outputs.

This is **not** a complete reproduction package for every manuscript result. It does not include the GSE191275 or GSE78523 source matrices or a current-five reanalysis of either cohort; the HPA source/version chain; the original cBioPortal response and complete TCGA mapping inputs; the historical 180-transcript discovery-rho generator; or the original GSE134520 QC object needed to reconcile the published 32,332-cell count. The reported limits remain in the manuscript.

## Human-data minimization

No raw sequencing files, cell-level matrices, participant-level clinical records, donor identifiers, ages, sexes, or sample-by-sample outputs are included in this release snapshot. The GSE163416 processed source tables are not redistributed; the manifest records their public GEO URLs, local byte counts, and SHA-256 checksums. Users should obtain source files directly from GEO and check applicable source terms.

The excluded local candidate files include GSE134520 donor metadata and its dependent donor summaries, as well as sample-level localization outputs. Their exclusion is intentional and does not change the manuscript's reported analysis or resolve the GSE134520 QC discrepancy.

## Files

- `scripts/rerun_gse163416_gene_level_20261009.R`: gene-level refit from three GSE163416 processed-count contrast tables.
- `scripts/gene_level_candidate_leave_one_sample_out_20261010.R`: leave-one-sample-out sensitivity analysis for the five candidates.
- `scripts/gene_level_go_enrichment_20261009.R`: GO-BP sensitivity enrichment.
- `results/Additional_file_2_validated_stage_monotonic_core_genes.csv`: historical comparator input.
- `results/Additional_file_9_gene_level_edger_results.csv`, `Additional_file_10_gene_level_sensitivity_candidates_cross_cohort.csv`, and `Additional_file_12_gene_level_GO_BP_sensitivity.csv`: aggregate results corresponding to manuscript supplements.
- `provenance/GSE163416_processed_input_manifest.json`: public source URLs and checksums for the processed inputs; input data are not bundled.
- `renv.lock` and `ENVIRONMENT_R-4.6.1_20261011.txt`: recorded software environment and package versions.

## Reproduce the included modules

Use R 4.6.1 and restore the recorded dependencies with `renv` from the repository root. Obtain the three GEO files named in the manifest and place them in one local input directory. The input filenames expected by the scripts are:

```text
GSE163416_N-CAG+IM.all.mRNA.txt.gz
GSE163416_N-Dys.all.mRNA.txt.gz
GSE163416_N-GC.all.mRNA.txt.gz
```

Run the primary refit:

```text
Rscript scripts/rerun_gse163416_gene_level_20261009.R <GSE163416_input_dir> results/Additional_file_2_validated_stage_monotonic_core_genes.csv <output_dir>
```

Run the five-candidate leave-one-sample-out analysis:

```text
Rscript scripts/gene_level_candidate_leave_one_sample_out_20261010.R <GSE163416_input_dir> <primary_output_dir>/gse163416_gene_level_edger_results.csv results/Additional_file_10_gene_level_sensitivity_candidates_cross_cohort.csv <output_dir>
```

Run the GO-BP analysis:

```text
Rscript scripts/gene_level_go_enrichment_20261009.R <primary_output_dir>/gse163416_gene_level_edger_results.csv results/Additional_file_10_gene_level_sensitivity_candidates_cross_cohort.csv <output_dir>
```

The recorded local reruns matched manuscript Additional files 9, 23, and 12, respectively. This documents the tested outputs and environment; it does not establish a fresh remote package restore or end-to-end reproduction of analyses whose source inputs are not included.

## Availability and reuse

This snapshot is public for inspection. **No software licence is included**, so no general permission to reuse, modify, or redistribute the code is granted by the authors in this release. GitHub's terms permit viewing and forking public repositories; that is not a substitute for a software licence. The licence status must be resolved before describing the code as freely reusable.

No persistent archive DOI is assigned in this snapshot. Cite the manuscript and repository version/commit when available; a permanent identifier should be added if a corresponding archive is created.

