nextflow run vcf2DCPCA.nf \
	--vcffile real-data/76g_chr1/chr1.vcf.gz \
	--output_dir real-data/76g_chr1/results \
	-resume \
	-with-report real-data/76g_chr1/results/`date +%Y%m%d_%H%M%S`_report.html \
	-with-dag real-data/76g_chr1/results/`date +%Y%m%d_%H%M%S`.DAG.html
