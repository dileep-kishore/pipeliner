#!/usr/bin/env nextflow

// BamStats for the bam file
process bam_stats{
    tag "$bam_rseqc_bamstats"

	executor = 'sge'
        clusterOptions = "-P ${params.project} -l h_rt=96:00:00 -l mem_total=5G -pe omp 12"

    publishDir "${params.outdir}/bam_stats", mode: 'copy'

    input:
    file bam_rseqc_bamstats

    output:
    file '*info.txt' into bam_stats_results

    script:
    """
    module load python
    module load rseqc/2.6.4
    bam_stat.py -i $bam_rseqc_bamstats > bam_stats_info.txt
    """
}

// GeneBodyCoverage for the reference (house-keeping) genes
process gene_body_coverage{
    tag "$bam_rseqc_genecoverage"

	executor = 'sge'
        clusterOptions = "-P ${params.project} -l h_rt=96:00:00 -l mem_total=5G -pe omp 12"

    publishDir "${params.outdir}/gene_coverage", mode: 'copy'

    input:
    file bam_rseqc_genecoverage
    file gtf from ref_gene_model

    output:
    file 'gene_coverage*' into gene_coverage_results
    stdout into gene_coverage_log

    // gtf must be the reference gene model
    // index file for bam must be available

    script:
    """
    module load python
    module load rseqc/2.6.4
    geneBody_coverage.py -r $gtf -i $bam_rseqc_genecoverage -o gene_coverage
    """
}

process junction_annotation{
    tag "$bam_rseqc_junc_annot"

	executor = 'sge'
        clusterOptions = "-P ${params.project} -l h_rt=96:00:00 -l mem_total=5G -pe omp 12"

    publishDir "${params.outdir}/junction_annotation", mode: 'copy'

    input:
    file bam_rseqc_junc_annot
    file gtf from ref_gene_model

    output:
    file 'junc_annot*' into junction_annotation_results

    // gtf must be the reference gene model
    // index file for bam must be available

    script:
    """
    module load python
    module load rseqc/2.6.4
    junction_annotation.py -i $bam_rseqc_junc_annot -o junc_annot -r $gtf
    """
}

// Mulitqc process that takes in all the output (hardcoded - needs to be changed)
process multiqc {
	executor = 'sge'
        clusterOptions = "-P ${params.project} -l h_rt=96:00:00 -l mem_total=5G -pe omp 12"

	publishDir "${params.outdir}/MultiQC", mode: 'copy'

    input:
    file ('fastqc/*') from fastqc_results.flatten().toList()
    file ('trimgalore/*') from trimgalore_results.flatten().toList()
    file ('alignment/*') from alignment_logs.flatten().toList()
    file ('stringtie/*') from stringtie_log.flatten().toList()
    file ('bam_stats/*') from bam_stats_results.flatten().toList()
    file ('gene_coverage/*') from gene_coverage_results.flatten().toList()
    file ('junction_annotation/*') from junction_annotation_results.flatten().toList()


    output:
    file "*multiqc_report.html"
    file "*multiqc_data"

    script:
    """
	 module load python/2.7.11
        module load multiqc/0.8

   multiqc -f ${params.outdir}/fastqc/ ${params.outdir}/STAR/ ${params.outdir}/trim_galore/ ${params.outdir}/stringtieFPKM/ ${params.outdir}/bam_stats/ ${params.outdir}/gene_coverage/ ${params.outdir}/junction_annotation/
    """
}
