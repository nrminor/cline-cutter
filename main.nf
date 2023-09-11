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
        .from( params.random_seeds )
        .splitCsv( header: false )
        .flatten()
    
    ch_sample_meta = Channel
        .fromPath( params.samplesheet )
	
	
	// Workflow steps 
	if ( params.precalled_vcf == "" ){

		DEMULTIPLEX_READS (
			ch_reads
		)

		INDEX_FOR_MAPPING ()

		MAP_WITH_BWA (
			DEMULTIPLEX_READS.out.flatten(),
			INDEX_FOR_MAPPING.out
		)

		CONVERT_AND_INDEX (
			MAP_WITH_BWA.out
		)

		VARIANT_CALL (
			CONVERT_AND_INDEX.out.bam.collect(),
			CONVERT_AND_INDEX.out.bai.collect()
		)

		RUN_DOWNSAMPLING (
			VARIANT_CALL.out,
			ch_sample_meta
		)

	} else {

		ch_vcf = Channel
			.fromPath( params.precalled_vcf )
		
		RUN_DOWNSAMPLING (
			ch_vcf,
			ch_sample_meta
		)

	}

    VCF_FILTERING (
        RUN_DOWNSAMPLING.out.vcf.flatten()
    )

    SNP_THINNING (
        VCF_FILTERING.out
    )

    FILTER_INDIVS (
        SNP_THINNING.out
    )

	CREATE_Q_PRIORS (
		FILTER_INDIVS.out
	)

    CONVERT_TO_MPGL (
        CREATE_Q_PRIORS.out
    )

    RUN_ENTROPY (
        ch_seeds,
        CONVERT_TO_MPGL.out
    )

    FIT_CLINE_MODELS (
        RUN_ENTROPY.out.groupTuple(),
        ch_sample_meta,
		RUN_DOWNSAMPLING.out.txt.collect()
    )
	
	
}
// --------------------------------------------------------------- //



// DERIVATIVE PARAMETER SPECIFICATION
// --------------------------------------------------------------- //
// Additional parameters that are derived from parameters set in nextflow.config

// --------------------------------------------------------------- //




// PROCESS SPECIFICATION 
// --------------------------------------------------------------- //

process DEMULTIPLEX_READS {
	
	/*
    This process does something described here
    */
	
	tag "${tag}"
	publishDir params.results, mode: 'copy'
	
	input:
	
	output:
    path "*.fastq.gz"
	
	script:
	"""
	"""
}

process INDEX_FOR_MAPPING {
	
	/*
    This process does something described here
    */
	
	tag "${tag}"
	publishDir params.results, mode: 'copy'
	
	input:
	
	output:
	
	script:
	"""
	"""
}

process MAP_WITH_BWA {
	
	/*
    This process does something described here
    */
	
	tag "${tag}"
	publishDir params.results, mode: 'copy'

    cpus 8
	
	input:
	
	output:
	
	script:
	"""
	"""
}

process CONVERT_AND_INDEX {
	
	/*
    This process does something described here
    */
	
	tag "${tag}"
	publishDir params.results, mode: 'copy'

    cpus 4
	
	input:
	
	output:
	
	script:
	"""
	"""
}

process VARIANT_CALL {
	
	/*
    This process does something described here
    */
	
	tag "${tag}"
	publishDir params.results, mode: 'copy'

    cpus 8
	
	input:
	
	output:
	
	script:
	"""
	"""
}

process RUN_DOWNSAMPLING {
	
	/*
    This process does something described here
    */
	
	tag "${tag}"
	publishDir params.results, mode: 'copy'

	cpus 4
	
	input:
	path vcf
	path samplesheet
	
	output:
	path "*.vcf", emit: vcf
	path "*.txt", emit: txt
	
	script:
	"""
	sample-by-coordinate.py \
	--vcf ${vcf} \
	--metadata ${samplesheet} \
	--distance_threshold 100 \
	--cores ${task.cpus} \
	--seed 14
	"""

}

process VCF_FILTERING {
	
	/*
    This process does something described here
    */
	
	tag "${tag}"
	publishDir params.results, mode: 'copy'
	
	input:
	path vcf
	
	output:
	path "*.vcf"
	
	script:
	"""
	vcftools --vcf ${vcf} --min-alleles 2 --max-alleles 2 \
    --maf 0.05 --max-missing 0.7 --recode --recode-INFO-all \
    --out ${params.project_name}
	"""

}

process SNP_THINNING {
	
	/*
    This process does something described here
    */
	
	tag "${tag}"
	publishDir params.results, mode: 'copy'
	
	input:
	path vcf
	
	output:
	path "*.vcf"
	
	script:
	"""
	vcftools --vcf ${vcf} --thin 20000 --recode --recode-INFO-all \
    --out ${params.project_name}_thinned
	"""

}

process FILTER_INDIVS {
	
	/*
    This process does something described here
    */
	
	tag "${tag}"
	publishDir params.results, mode: 'copy'
	
	input:
	path vcf
	
	output:
	path "*.vcf"
	
	script:
	"""
	vcftools --vcf ${vcf} --missing-indv
    awk '$5 >= 0.7 {print $1}' out.imiss > individuals_to_remove.txt
    vcftools --vcf ${vcf} --remove individuals_to_remove.txt \
	--recode --recode-INFO-all \
    --out ${params.project_name}_final
	"""

}

process CREATE_Q_PRIORS {
	
	/*
    This process does something described here
    */
	
	tag "${tag}"
	publishDir params.results, mode: 'copy'

    cpus 1
	
	input:
	path vcf
	
	output:
	tuple path(vcf), path("*.txt")
	
	shell:
	'''
	N=${bcftools query -l !{vcf} | wc -l}
	touch starting_q.txt
	for (( i=1; i<=$num_rows; i++ ))
	do
		echo "0.5" >> "starting_q.txt"
	done
	'''

}

process CONVERT_TO_MPGL {
	
	/*
    This process does something described here
    */
	
	tag "${tag}"
	publishDir params.results, mode: 'copy'

    cpus 1
	
	input:
	tuple path(vcf), path(starting_q)
	
	output:
	tuple path("*.mpgl"), path(starting_q)
	
	script:
	"""
	vcf2mpgl.R ${vcf}
	"""

}

process RUN_ENTROPY {
	
	/*
    This process does something described here
    */
	
	tag "${tag}"
	publishDir params.results, mode: 'copy'

    cpus 1
	
	input:
    each val(random_seed)
	tuple path(mpgl), path(starting_q)
	
	output:
	tuple val(subsample), path("*.hdf5")
	
	script:
	subsample = file(mpgl.toString()).getSimpleName()
	"""
	entropy -i ${mpgl} \
    -r ${random_seed} -q ${starting_q} \
    -m 1 -n 2 -k 2 -w 1 -Q 1 -l 120000 -b 30000 -t 30 \
    -o ${params.project_name}_${seed}.hdf5
	"""

}

process FIT_CLINE_MODELS {
	
	/*
    This process does something described here
    */
	
	tag "${tag}"
	publishDir params.results, mode: 'copy'

    cpus 1
	
	input:
	tuple val(subsample), tuple(path(hdf5_1), path(hdf5_2), path(hdf5_3)) 
    path samplesheet
	path subset_files
	
	output:
	path "*.pdf"
	
	script:
	subset_file = "${subsample}_sample.txt"
	"""
	cline_fitting.R ${hdf5}
	"""

}

// --------------------------------------------------------------- //