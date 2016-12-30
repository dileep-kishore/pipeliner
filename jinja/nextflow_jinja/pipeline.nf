#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
RNA-Seq Analysis Pipeline
*/

params.reads = "$baseDir/data/ggal/*_{1,2}.fq"
params.annot = "$baseDir/data/ggal/ggal_1_48850000_49020000.bed.gff"
params.genome = "$baseDir/data/ggal/ggal_1_48850000_49020000.Ggal71.500bpflank.fa"
params.outdir = "results"

log.info "R N A T O Y   P I P E L I N E    "
log.info "================================="
log.info "genome             : ${params.genome}"
log.info "annotat            : ${params.annot}"
log.info "reads              : ${params.reads}"
log.info "outdir             : ${params.outdir}"

/*
 * the reference genome file
 */
genome_file = file(params.genome)
annotation_file = file(params.annot)

/*
 * Create the `read_pairs` channel that emits tuples containing three elements:
 * the pair ID, the first read-pair file and the second read-pair file
 */
Channel
    .fromFilePairs( params.reads )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .set { read_pairs }

workflow.onComplete {
	println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}

/*
 * Step 1. Builds the genome index required by the mapping process
 */
process buildIndex {
    tag "$genome_file.baseName"

    input:
    file genome_file

    output:
    file 'genome.index*' into genome_index

    """
    bowtie2-build --threads ${task.cpus} ${genome_file} genome.index
    """
}

/*
 * Step 2. Maps each read-pair by using Tophat2 mapper tool
 */
process mapping {
    tag "$pair_id"

    input:
    file 'genome.index.fa' from genome_file
    file annotation_file
    file genome_index from genome_index.first()
    set pair_id, file(reads) from read_pairs

    output:
    set pair_id, "accepted_hits.bam" into bam

    """
    tophat2 -p ${task.cpus} --GTF $annotation_file genome.index ${reads}
    mv tophat_out/accepted_hits.bam .
    """
}

/*
 * Step 3. Assembles the transcript by using the "cufflinks" tool
 */
process makeTranscript {
    tag "$pair_id"
    publishDir params.outdir, mode: 'copy'

    input:
    file 'anno.gtf' from annotation_file
    set pair_id, file(bam_file) from bam

    output:
    set pair_id, file('transcript_*.gtf') into transcripts

    """
    cufflinks --no-update-check -q -p ${task.cpus} -G anno.gtf ${bam_file}
    mv transcripts.gtf transcript_${pair_id}.gtf
    """
}

