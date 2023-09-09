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

    RUN_ADMIXTURE (
        RUN_DOWNSAMPLING.out.flatten()
    )

    VCF_FILTERING (
        RUN_DOWNSAMPLING.out.flatten()
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
        FILTER_INDIVS.out
    )

    RUN_ENTROPY (
        ch_seeds,
        CONVERT_TO_MPGL.out,
        CREATE_Q_PRIORS.out.collect()
    )

    FIT_CLINE_MODELS (
        RUN_ENTROPY.out,
        ch_sample_meta
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
	path "*.vcf"
	
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

process RUN_ADMIXTURE {
	
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
    vcftools --vcf $vcf --remove individuals_to_remove.txt --recode --recode-INFO-all \
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
	path "*.mpgl"
	
	shell:
	'''
	N=${bcftools query -l !{vcf} | wc -l}
	touch starting_q.txt
	for (( i=1; i<=$num_rows; i++ ))
	do
		echo "0.5" >> "$output_file"
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
	path vcf
	
	output:
	path "*.mpgl"
	
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
	path mpgl
    path starting_qs
	
	output:
	path "*.hdf5"
	
	script:
	subsample = file(mpgl.toString()).getSimpleName()
	q_file = "${subsample}_q.txt"
	"""
	entropy -i ${mpgl} \
    -r ${random_seed} -q ${q_file} \
    -m 1 -n 2 -k 2 -w 1 -Q 1 -l 120000 -b 30000 -t 30 \
    -o ${params.project_name}.hdf5
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
	path hdf5
    path samplesheet
	
	output:
	path "*.pdf"
	
	script:
	"""
	cline_fitting.R ${hdf5}
	"""

}

// --------------------------------------------------------------- //