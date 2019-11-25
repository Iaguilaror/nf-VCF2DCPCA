#!/usr/bin/env nextflow

/*================================================================
The MORETT LAB presents...

  The VCF to PCA pipeline

-

==================================================================
Version: 0.0.1
Project repository: PENDING
==================================================================
Authors:

- Bioinformatics Design
 Israel Aguilar-Ordonez (iaguilaror@gmail)

- Bioinformatics Development
 Israel Aguilar-Ordonez (iaguilaror@gmail)

- Nextflow Port
 Israel Aguilar-Ordonez (iaguilaror@gmail)

=============================
Pipeline Processes In Brief:

Pre-processing:

Core-processing:
  _001_filterSNP
	_002_recodeGTdata
	_003_DCPCA

Pos-processing

================================================================*/

/* Define the help message as a function to call when needed *//////////////////////////////
def helpMessage() {
	log.info"""
  ==========================================
  PENDING
  - PENDING
  v${version}
  ==========================================

	Usage:

  nextflow run vcf2DCPCA.nf --vcffile <path to input 1> [--output_dir path to results ]

	  --vcffile    <- compressed vcf file for annotation;
				vcf file must be annotated with https://github.com/Iaguilaror/nf-VEPextended
				accepted extension is vcf.gz;
	  --output_dir     <- directory where results, intermediate and log files will bestored;
				default: same dir where --query_fasta resides
	  -resume	   <- Use cached results if the executed project has been run before;
				default: not activated
				This native NF option checks if anything has changed from a previous pipeline execution.
				Then, it resumes the run from the last successful stage.
				i.e. If for some reason your previous run got interrupted,
				running the -resume option will take it from the last successful pipeline stage
				instead of starting over
				Read more here: https://www.nextflow.io/docs/latest/getstarted.html#getstart-resume
	  --help           <- Shows Pipeline Information
	  --version        <- Show pipeline version
	""".stripIndent()
}

/*//////////////////////////////
  Define pipeline version
  If you bump the number, remember to bump it in the header description at the begining of this script too
*/
version = "0.0.1"

/*//////////////////////////////
  Define pipeline Name
  This will be used as a name to include in the results and intermediates directory names
*/
pipeline_name = "VCF2DCPCA"

/*
  Initiate default values for parameters
  to avoid "WARN: Access to undefined parameter" messages
*/
params.vcffile = false  //if no inputh path is provided, value is false to provoke the error during the parameter validation block
params.help = false //default is false to not trigger help message automatically at every run
params.version = false //default is false to not trigger version message automatically at every run

/*//////////////////////////////
  If the user inputs the --help flag
  print the help message and exit pipeline
*/
if (params.help){
	helpMessage()
	exit 0
}

/*//////////////////////////////
  If the user inputs the --version flag
  print the pipeline version
*/
if (params.version){
	println "VCF 2 DCPCA v${version}"
	exit 0
}

/*//////////////////////////////
  Define the Nextflow version under which this pipeline was developed or successfuly tested
  Updated by iaguilar at FEB 2019
*/
nextflow_required_version = '18.10.1'
/*
  Try Catch to verify compatible Nextflow version
  If user Nextflow version is lower than the required version pipeline will continue
  but a message is printed to tell the user maybe it's a good idea to update her/his Nextflow
*/
try {
	if( ! nextflow.version.matches(">= $nextflow_required_version") ){
		throw GroovyException('Your Nextflow version is older than Pipeline required version')
	}
} catch (all) {
	log.error "-----\n" +
			"  This pipeline requires Nextflow version: $nextflow_required_version \n" +
      "  But you are running version: $workflow.nextflow.version \n" +
			"  The pipeline will continue but some things may not work as intended\n" +
			"  You may want to run `nextflow self-update` to update Nextflow\n" +
			"============================================================"
}

/*//////////////////////////////
  INPUT PARAMETER VALIDATION BLOCK
/* Check if inpuths were provided
 * if they were not provided, they keep the 'false' value assigned in the parameter initiation block above and this test fails
*/
if ( !file(params.vcffile).exists() ) {
  log.error "Input file does not exist\n\n" +
  "Please provide valid path for the --vcffile and make sure that it has an index \n\n" +
  "For more information, execute: nextflow run vcf2DCPCA.nf --help"
  exit 1
}

/*  Check that extension of input 1 is .vcf.gz
 * to understand regexp use, and '==~' see https://www.nextflow.io/docs/latest/script.html#regular-expressions
*/
if ( !(file(params.vcffile).getName() ==~ /.+\.vcf\.gz$/) ) {
	log.error " --vcffile must have .vcf.gz extension \n\n" +
	"For more information, execute: nextflow run vcf2DCPCA.nf --help"
  exit 1
}

/*
Output directory definition
Default value to create directory is the parent dir of --vcffile
*/
params.output_dir = file(params.vcffile).getParent()

/*
  Results and Intermediate directory definition
  They are always relative to the base Output Directory
  and they always include the pipeline name in the variable (pipeline_name) defined by this Script

  This directories will be automatically created by the pipeline to store files during the run
*/
results_dir = "${params.output_dir}/${pipeline_name}-results/"
intermediates_dir = "${params.output_dir}/${pipeline_name}-intermediate/"

/*//////////////////////////////
  LOG RUN INFORMATION
*/
log.info"""
==========================================
The VCF 2 DCPCA
- PENDING
v${version}
==========================================
"""
log.info "--Nextflow metadata--"
/* define function to store nextflow metadata summary info */
def nfsummary = [:]
/* log parameter values beign used into summary */
/* For the following runtime metadata origins, see https://www.nextflow.io/docs/latest/metadata.html */
nfsummary['Resumed run?'] = workflow.resume
nfsummary['Run Name']			= workflow.runName
nfsummary['Current user']		= workflow.userName
/* string transform the time and date of run start; remove : chars and replace spaces by underscores */
nfsummary['Start time']			= workflow.start.toString().replace(":", "").replace(" ", "_")
nfsummary['Script dir']		 = workflow.projectDir
nfsummary['Working dir']		 = workflow.workDir
nfsummary['Current dir']		= workflow.launchDir
nfsummary['Launch command'] = workflow.commandLine
log.info nfsummary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "\n\n--Pipeline Parameters--"
/* define function to store nextflow metadata summary info */
def pipelinesummary = [:]
/* log parameter values beign used into summary */
pipelinesummary['VCFfile']			= params.vcffile
pipelinesummary['Results Dir']		= results_dir
pipelinesummary['Intermediate Dir']		= intermediates_dir
/* print stored summary info */
log.info pipelinesummary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "==========================================\nPipeline Start"

/*//////////////////////////////
  PIPELINE START
*/

/*
	READ INPUTS
*/

/* Load vcf files AND TABIX INDEX into channel */
Channel
  .fromPath("${params.vcffile}*")
	.toList()
  .set{ vcf_inputs }

/* _001_filterSNP */
/* Read mkfile module files */
Channel
	.fromPath("${workflow.projectDir}/mkmodules/mk-filterSNP/*")
	.toList()
	.set{ mkfiles_001 }

process _001_filterSNP {
	echo true

	publishDir "${intermediates_dir}/_001_filterSNP/",mode:"symlink"

	input:
	file sample from vcf_inputs
	file mk_files from mkfiles_001

	output:
	file "*.vcf" into results_001, also_results_001

	"""
	bash runmk.sh
	"""

}

/* _002_recodeGTdata */
/* Read mkfile module files */
Channel
	.fromPath("${workflow.projectDir}/mkmodules/mk-recodeGTdata/*")
	.toList()
	.set{ mkfiles_002 }

process _002_recodeGTdata {
	echo true

	publishDir "${intermediates_dir}/_002_recodeGTdata/",mode:"symlink"

	input:
	file sample from results_001
	file mk_files from mkfiles_002

	output:
	file "*.txt" into results_002

	"""
	bash runmk.sh
	"""

}

/* _003_DCPCA */
/* Read mkfile module files */
Channel
	.fromPath("${workflow.projectDir}/mkmodules/mk-DCPCA/*")
	.toList()
	.set{ mkfiles_003 }

process _003_DCPCA {
	// echo true

	publishDir "${results_dir}/_003_DCPCA/",mode:"copy"

	input:
	file sample from results_002
	file vcf from also_results_001
	file mk_files from mkfiles_003

	output:
	file "*" into results_003

	"""
	bash runmk.sh
	"""

}
