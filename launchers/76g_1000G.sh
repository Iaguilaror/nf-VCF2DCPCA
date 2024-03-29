nextflow run vcf2DCPCA.nf \
	--vcffile real-data/76g_1000genomes/76g_1000GP-population_set.vcf.gz \
	--output_dir real-data/76g_1000genomes/results \
	-resume \
	-with-report real-data/76g_1000genomes/results/`date +%Y%m%d_%H%M%S`_report.html \
	-with-dag real-data/76g_1000genomes/results/`date +%Y%m%d_%H%M%S`.DAG.html
