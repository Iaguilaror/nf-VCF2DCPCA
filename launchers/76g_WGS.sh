nextflow run vcf2DCPCA.nf \
	--vcffile real-data/76g_autosomes/76g_PASS_ANeqorgt150_autosomes_and_XY.filtered.vcf.gz \
	--output_dir real-data/76g_autosomes/results \
	-resume \
	-with-report real-data/76g_autosomes/results/`date +%Y%m%d_%H%M%S`_report.html \
	-with-dag real-data/76g_autosomes/results/`date +%Y%m%d_%H%M%S`.DAG.html
