#!/usr/bin/env nextflow

nextflow.enable.dsl = 2



// WORKFLOW SPECIFICATION
// --------------------------------------------------------------- //
workflow {

	// input channels
    ch_reads = Channel
        .fromPath( "${params.input_dir}/*.fastq.gz" )
		.collect()

    ch_seeds = Channel
        .of( 1, 2, 3 )

	// ch_proportions = Channel
	// 	.of( 0.5, 0.8, 0.9 )

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

	// ch_sample_meta
	// 	.combine( ch_proportions )
	// 	.view()

	RUN_DOWNSAMPLING (
		FILTER_INDIVS.out,
		ch_sample_meta
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
        ch_seeds,
        CREATE_Q_PRIORS.out
    )

    FIT_CLINE_MODELS (
        RUN_ENTROPY.out
			.groupTuple( sort: true )
			.map { sample, hdf5s -> tuple( sample, hdf5s[0], hdf5s[1], hdf5s[2] ) },
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
	"""
	vcftools --vcf ${vcf} --min-alleles 2 --max-alleles 2 \
    --maf 0.05 --max-missing 0.7 --recode --recode-INFO-all \
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
	vcftools --vcf ${vcf} --thin 20000 --recode --recode-INFO-all \
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

	publishDir params.downsampled, mode: 'copy'

	cpus 4
	time '6h'

	input:
	path vcf
	path samplesheet

	output:
	path "*.vcf", emit: vcf

	script:
	"""
	sample-by-coordinate.py \
	--vcf ${vcf} \
	--proportion ${params.proportion} \
	--metadata ${samplesheet} \
	--distance_threshold ${params.distance_threshold} \
	--cores ${task.cpus} \
	--seed 14
	"""

}

process RECORD_FINAL_ROSTER {

	/*
	*/

	publishDir params.downsampled, mode: 'copy'

	cpus 1
	time '10m'

	input:
	path vcf

	output:
	path "*.txt"

	script:
	downsampling_regime = file(vcf.toString()).getSimpleName()
	"""
	bcftools query -l ${vcf} > ${downsampling_regime}.txt
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

	tag "${sample_regime}"

    cpus 1
	time '10m'

	input:
	path mpgl

	output:
	tuple path(mpgl), path("*.txt")

	script:
	sample_regime = file(mpgl.toString()).getSimpleName().replace(".recode", "")
	"""
	mpgl_sample_size.py -m ${mpgl} -q ${params.starting_q} -l ${sample_regime}
	"""

}

process RUN_ENTROPY {

	/* */

	tag "${subsample}, ${random_seed}"
	publishDir params.entropy, mode: 'copy'

    cpus 8
	time '7d'

	input:
    each random_seed
	tuple path(mpgl), path(starting_q)

	output:
	tuple val(subsample), path("*.hdf5")

	script:
	subsample = file(mpgl.toString()).getSimpleName()
	"""
	entropy -i ${mpgl} \
    -r ${random_seed} -q ${starting_q} \
    -m 1 -n 2 -k 2 -w 1 -Q 1 -l 120000 -b 30000 -t 30 \
    -o ${subsample}_${random_seed}.hdf5
	"""

}

// process VISUALIZE_ENTROPY_TRACE {}

process FIT_CLINE_MODELS {

	/* */

	tag "${subsample}"
	publishDir params.clines, mode: 'copy'

    cpus 1
	time '8h'

	input:
	tuple val(subsample), path(hdf5_1), path(hdf5_2), path(hdf5_3)
    path metadata_files

	output:
	path "*", emit: all_files
	tuple val(subsample), path("${subsample}_model_logs.txt"), path("*_aic.tsv"), emit: modeling_logs

	script:
	"""
	cline_fitting.R &> ${subsample}_model_logs.txt
	"""

}

process EVAL_MODEL_PERFORMANCE {

	/* */

	tag "${subsample}"
	publishDir params.clines, mode: 'copy'

    cpus 1
	time '10m'

	input:
	tuple val(subsample), path(modeling_logs), path(aic_values)

	output:
	path "*"

	script:
	"""
	collate_model_evals.py "${subsample}" "${subsample}_model_logs.txt"
	"""


}

// process REPORT_MODEL_PERFORMANCE {}

// --------------------------------------------------------------- //
