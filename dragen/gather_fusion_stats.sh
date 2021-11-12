#!/bin/bash


find ./samples -type f | ~/code/FusionBenchmarking/util/make_file_listing_input_table.pl > fusion_result_file_listing.dat

~/code/FusionBenchmarking/benchmarking/collect_preds.pl fusion_result_file_listing.dat > preds.collected

~/code/FusionBenchmarking/benchmarking/map_gene_symbols_to_gencode.pl preds.collected ~/code/FusionBenchmarking/resources/genes.coords.gz ~/code/FusionBenchmarking/resources/genes.aliases > preds.collected.gencode_mapped 2> preds.collected.gencode_mapped.warnings

../STAR-Fusion.v1.10.1/FusionAnnotator/FusionAnnotator --genome_lib_dir ../GRCh37_gencode_v19_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir/ --annotate preds.collected.gencode_mapped -C 2 > preds.collected.gencode_mapped.wAnnot

~/code/FusionBenchmarking/benchmarking/filter_collected_preds.pl preds.collected.gencode_mapped.wAnnot > preds.collected.gencode_mapped.wAnnot.filt

~/code/FusionBenchmarking/benchmarking/fusion_preds_to_TP_FP_FN.pl --truth_fusions ~/code/FusionBenchmarking/simulated_data/sim_101/sim_101.truth_set.dat --fusion_preds preds.collected.gencode_mapped.wAnnot.filt --allow_reverse_fusion --allow_paralogs ~/code/FusionBenchmarking/resources/paralog_clusters.dat > preds.collected.gencode_mapped.wAnnot.filt.scored

~/code/FusionBenchmarking/benchmarking/all_TP_FP_FN_to_ROC.pl preds.collected.gencode_mapped.wAnnot.filt.scored > preds.collected.gencode_mapped.wAnnot.filt.scored.ROC

~/code/FusionBenchmarking/benchmarking/calc_PR.py --in_ROC preds.collected.gencode_mapped.wAnnot.filt.scored.ROC --out_PR preds.collected.gencode_mapped.wAnnot.filt.scored.PR.AUC| sort -k2,2gr | tee preds.collected.gencode_mapped.wAnnot.scored.PR.AUC
