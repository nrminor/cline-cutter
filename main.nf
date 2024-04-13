#!/usr/bin/env nextflow

nextflow.enable.dsl = 2



// WORKFLOW SPECIFICATION
// --------------------------------------------------------------- //

// prints to the screen and to the log
// shout out to the following for the ascii art:
// https://patorjk.com/software/taag/#p=testall&f=Graffiti&t=Cline%20Cutter
log.info	"""

			..######..##.......####.##....##.########.....######..##.....##.########.########.########.########.
			.##....##.##........##..###...##.##..........##....##.##.....##....##.......##....##.......##.....##
			.##.......##........##..####..##.##..........##.......##.....##....##.......##....##.......##.....##
			.##.......##........##..##.##.##.######......##.......##.....##....##.......##....######...########.
			.##.......##........##..##..####.##..........##.......##.....##....##.......##....##.......##...##..
			.##....##.##........##..##...###.##..........##....##.##.....##....##.......##....##.......##....##.
			..######..########.####.##....##.########.....######...#######.....##.......##....########.##.....##

			CLINE-CUTTER:
			A Master's Project from University of Wyoming Exploring How Space and Sampling
			Affect Cline Modeling across an Avian Hybrid Zone
			(version 0.1.0)
			===============================================================================

			Inputs:
			------------------------------------------
			Project name       : ${params.project_name}
			Samplesheet        : ${params.samplesheet}
			Pre-called VCF     : ${params.precalled_vcf}

			Variant Filtering Settings:
			-------------------------------------------
			Minor allele freq. : ${params.maf}
			Max missingness    : ${params.max_missing}
			HWE cutoff         : ${params.hwe}
			LD thinning cutoff : ${params.thinning}

			Downsampling Settings:
			-------------------------------------------
			Proportions        : ${params.proportions}
			Seeds              : ${params.seeds}
			Distance threshold : ${params.distance_threshold}

			Modeling Settings:
			-------------------------------------------
			Hybrid index prior : ${params.starting_q}
			Seeds              : ${params.seeds}

			"""
			.stripIndent()


workflow {

	assert params.precalled_vcf.endsWith(".vcf") : "Input VCF must be uncompressed and end with '.vcf'."

	// input channels
    ch_seeds = Channel
        .of( params.seeds )
		.splitCsv( header: false )
		.flatten()
		.map{ x -> x.stripIndent().toInteger() }

	ch_proportions = Channel
        .of( params.proportions )
		.splitCsv( header: false )
		.flatten()
		.map{ x -> x.stripIndent().toFloat() }

    ch_sample_meta = Channel
        .fromPath( params.samplesheet )

	ch_vcf = Channel
		.fromPath( params.precalled_vcf )


	// Processes
	VCF_FILTERING (
		ch_vcf
	)

    SNP_THINNING (
        VCF_FILTERING.out
    )

    FILTER_INDIVS (
        SNP_THINNING.out
    )

	RUN_DOWNSAMPLING (
		FILTER_INDIVS.out,
		ch_sample_meta
			.combine(
				ch_proportions
					.combine( ch_seeds )
			)
	)

	RECORD_FINAL_ROSTER (
		RUN_DOWNSAMPLING.out.vcf.flatten()
	)

    CONVERT_TO_MPGL (
        RUN_DOWNSAMPLING.out.vcf.collect()
    )

	CREATE_Q_PRIORS (
		CONVERT_TO_MPGL.out.flatten()
	)

    RUN_ENTROPY (
        CREATE_Q_PRIORS.out
			.combine(ch_seeds)
    )

    FIT_CLINE_MODELS (
        RUN_ENTROPY.out
			.groupTuple( by: [0,1,2], sort: true )
			.map { sample, seed, hdf5s -> tuple( sample, seed, hdf5s[0], hdf5s[1], hdf5s[2] ) },
        ch_sample_meta
			.mix(RECORD_FINAL_ROSTER.out.collect())
			.collect()
    )

	EVAL_MODEL_PERFORMANCE (
		FIT_CLINE_MODELS.out.modeling_logs
	)


}
// --------------------------------------------------------------- //



// DERIVATIVE PARAMETER SPECIFICATION
// --------------------------------------------------------------- //
// Additional parameters that are derived from parameters set in nextflow.config

// VCF files
params.variants = params.results + "/1_VCF_files"

// VCFs sub-results
params.filtered = params.variants + "/1_filtered"
params.thinned = params.variants + "/2_thinned"
params.no_missing = params.variants + "/3_no_missing_indiv"
params.downsampled = params.variants + "/4_downsampled"

// Analyses and visualizations
params.analyses = params.results + "/2_analyses"
params.entropy = params.analyses + "/1_entropy"
params.clines = params.analyses + "/2_clines"

// --------------------------------------------------------------- //




// PROCESS SPECIFICATION
// --------------------------------------------------------------- //

process VCF_FILTERING {

	/* */

	publishDir params.filtered, mode: 'copy'

	time '6h'

	input:
	path vcf

	output:
	path "*.vcf"

	script:
	subsample = file(vcf.toString()).getSimpleName().replace("_sample", "")
	assert vcf.toString().endsWith(".vcf") : "Input VCF must be uncompressed and end with '.vcf'."
	"""
	vcftools \
	--vcf ${vcf} \
	--min-alleles 2 --max-alleles 2 \
    --maf ${params.maf} --max-missing ${params.max_missing} --hwe ${params.hwe} \
	--recode --recode-INFO-all \
    --out ${params.project_name}_${subsample}
	"""

}

process SNP_THINNING {

	/* */

	publishDir params.thinned, mode: 'copy'

	time '1h'

	input:
	path vcf

	output:
	tuple path("*.vcf"), val(simple_name)

	script:
	simple_name = file(vcf.toString()).getSimpleName()
	"""
	vcftools \
	--vcf ${vcf} \
	--thin ${params.thinning} \
	--recode --recode-INFO-all \
    --out ${simple_name}_thinned
	"""

}

process FILTER_INDIVS {

	/* */

	publishDir params.no_missing, mode: 'copy'

	time '1h'

	input:
	tuple path(vcf), val(simple_name)

	output:
	path "*.vcf"

	shell:
	'''
	# find individuals where all sites are marked as missing
	vcftools --vcf !{vcf} --missing-indv

	# parse those into a simple text file
    awk '$5 >= 0.7 {print $1}' out.imiss > individuals_to_remove.txt

	# use vcftools again to remove individuals in the all-missing txt file
	vcftools \
	--vcf !{vcf} \
	--remove individuals_to_remove.txt \
	--recode --recode-INFO-all \
    --out !{simple_name}_no_missing
	'''

}

process RUN_DOWNSAMPLING {

	/* */

	tag "${proportion}, ${seed}"
	publishDir params.downsampled, mode: 'copy'

	cpus 4
	time '6h'

	input:
	each path(vcf)
	tuple path(samplesheet), val(proportion), val(seed)

	output:
	path "*.vcf", emit: vcf

	script:
	"""
	sample-by-coordinate.py \
	--vcf ${vcf} \
	--proportion ${proportion} \
	--metadata ${samplesheet} \
	--distance_threshold ${params.distance_threshold} \
	--cores ${task.cpus} \
	--seed ${seed}
	"""

}

process RECORD_FINAL_ROSTER {

	/*
	*/

	publishDir params.downsampled, mode: 'copy'

	cpus 1
	time '10m'

	input:
	path(vcf)

	output:
	path "*.txt"

	script:
	file_base = file(vcf.toString()).getSimpleName().split(".recode")[0]
	name_parts = file_base.split("_")
	downsampling_regime = name_parts[0]
	proportion = name_parts[1]
	seed = name_parts[2]
	"""
	bcftools query -l ${vcf} > ${downsampling_regime}_${proportion}_${seed}.txt
	"""

}

process CONVERT_TO_MPGL {

	/* */

	publishDir params.entropy, mode: 'copy'

    cpus 12
	time '1h'

	input:
	path vcf_files

	output:
	path "*.mpgl"

	script:
	"""
	vcf2mpgl.jl
	"""

}

process CREATE_Q_PRIORS {

	/* */

	tag "${downsampling_regime}, ${proportion}, ${seed}"

    cpus 1
	time '10m'

	input:
	path mpgl

	output:
	tuple path(mpgl), path("*.txt"), val(proportion), val(seed)

	script:
	file_base = file(mpgl.toString()).getBaseName().split(".recode")[0]
	name_parts = file_base.replace("_sample", "").split("_")
	downsampling_regime = name_parts[0]
	proportion = name_parts[1]
	seed = name_parts[2]
	"""
	mpgl_sample_size.py \
	-m ${mpgl} \
	-q ${params.starting_q} \
	-l ${downsampling_regime}_${proportion}_${seed}
	"""

}

process RUN_ENTROPY {

	/* */

	tag "${subsample}, ${proportion}, ${seed}, ${entropy_seed}"
	publishDir "${params.entropy}/${subsample}", mode: 'copy'

    cpus 8
	time '7d'

	input:
	tuple path(mpgl), path(starting_q), val(proportion), val(seed), val(entropy_seed)

	output:
	tuple val(subsample), val(proportion), val(seed), path("${subsample}_${proportion}_${seed}_${entropy_seed}.hdf5")

	script:
	subsample = file(mpgl.toString()).getSimpleName().split("_")[0]
	"""
	entropy -i ${mpgl} \
    -r ${entropy_seed} -q ${starting_q} \
    -m 1 -n 2 -k 2 -w 1 -Q 1 -l 120000 -b 30000 -t 30 \
    -o ${subsample}_${proportion}_${seed}_${entropy_seed}.hdf5
	"""

}

// process VISUALIZE_ENTROPY_TRACE {}

process FIT_CLINE_MODELS {

	/* */

	tag "${subsample}, ${proportion}, ${seed}"
	publishDir "${params.clines}/${subsample}", mode: 'copy'

    cpus 1
	time '8h'

	input:
	tuple val(subsample), val(proportion), val(seed), path(hdf5_1), path(hdf5_2), path(hdf5_3)
    path metadata_files

	output:
	path "*", emit: all_files
	tuple val(subsample), val(proportion), val(seed), path("${subsample}_model_logs.txt"), path("*_aic.tsv"), emit: modeling_logs

	script:
	"""
	cline_fitting.R &> ${subsample}_${seed}_model_logs.txt
	"""

}

process EVAL_MODEL_PERFORMANCE {

	/* */

	tag "${subsample}, ${proportion}, ${seed}"
	publishDir "${params.clines}/${subsample}", mode: 'copy'

    cpus 1
	time '10m'

	input:
	tuple val(subsample), val(proportion), val(seed), path(modeling_logs), path(aic_values)

	output:
	path "*"

	script:
	"""
	collate_model_evals.py "${subsample}_${proportion}_${seed}" "${subsample}_model_logs.txt"
	"""


}

// process REPORT_MODEL_PERFORMANCE {}

// --------------------------------------------------------------- //
