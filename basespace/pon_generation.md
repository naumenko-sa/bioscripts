# How to generate a PON for EMG for CNV calling with Dragen for a panel assay

1. Investigate the variability of the coverage in assay - see cheo/pon_targets.qmd template.
2. Select samples with median_cov over genes > 300X, min_cov > 150X. The maximum N of samples in BS for a PON = 100.
3. Calculate PON with Basespace/ DRAGEN baseline Builder v4.3.8
3.1 Create project
3.2 upload bams
one command produces one dataset with multiple files
bs dataset upload -p [project_id] --type common.files --exclude '*' --include '*.bam' --dry-run --recursive .
note that sample name will be taken from the header @RG/SM tag no duplicates allowed in the dataset
3.3 upload bed
bs dataset upload -p [project_id] --type common.files PANELTargets_EMG_v2_hg38_2025-05-07.bed
3.3 bed file should be the same for BSSH and EMG - 4 column format, no header

APP settings (ran v4.3.8 - has to match Emedgene Dragen)
- select project
- Baseline Mode: CNV
- select bam files
- select reference hg38 Multigenome/Pangenome
- select bed file

EMG specific settings
- Switch ON: Ignore Duplicate Reads in CNV Baseline Files"
- change "Generated Combined Counts file for CNV" to "GC Corrected"

4. Download pon file
bs list appsessions  --project-name [PROJECT_NAME]
bs list datasets  --filter-field AppSession.Id --filter-term [session_id] | grep pon
bs download dataset -i [ds.id] -o .
pon is in pon.combined.counts.txt.gz

5. Plot correlations, select samples with Rsq >0.95.

6. Attach PON to EMG:
https://help.emg.illumina.com/release-notes/workbench-and-pipeline-updates/new_in_emedgene_37_-february_17_2025#attach-pon-to-kit

Male/Female ratio should be ~ 50/50
determine sex for example with https://github.com/adrianodemarino/Determine_sex_from_bam
(in targets we don't have probes on Y)
It also visible on correlation heatmap.

Cost for 89 Panel samples ~ 12 iCredits

If post-analysis of a PON reveals suboptimal samples, 
a better subset could be selected without re-running the PON generation algorithm with pon_subset.py script.
If the BED file remains the same, the coverage values don't differ between reruns of the same sample.