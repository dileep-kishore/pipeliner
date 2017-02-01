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

workflow.onComplete {
	println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}

process fastqc {
/* executor = 'sge'
   * clusterOptions = "-P ${PROJECT} -l h_rt=24:00:00 -l mem_total=5G"
     memory { 2.GB * task.attempt }
   * time { 4.h * task.attempt }
   * errorStrategy { task.exitStatus == 143 ? 'retry' : 'ignore' }
   * maxRetries 3
   * maxErrors '-1'
*/
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
    file(reads:'*') from read_files_fastqc

    output:
    file '*_fastqc.{zip,html}' into fastqc_results

    """
   /restricted/projectnb/pulmseq/kkarri_netflow_test/FastQC/fastqc -q $reads
    """
}

process trim_galore {

  /*  executor = 'sge'
*	clusterOptions = "-P ${PROJECT} -l h_rt=24:00:00 -l mem_total=5G"
*	cpus 3
 *   memory { 3.GB * task.attempt }
  *  time { 16.h * task.attempt }
   * errorStrategy { task.exitStatus == 143 ? 'retry' : 'terminate' }
    *maxRetries 3
   * maxErrors '-1'
*/
    publishDir "${params.outdir}/trim_galore", mode: 'copy'

    input:
    file(reads:'*') from read_files_trimming

    output:
    file '*fq.gz' into trimmed_reads
    file '*trimming_report.txt' into trimgalore_results

    script:
    single = reads instanceof Path
    if(single) {
	  """
	module load python2.7/Python-2.7.3_gnu446
	module load cutadapt/1.7.1_Python-2.7.3
     /restricted/projectnb/pulmseq/kkarri_netflow_test/trim_galore --gzip $reads
        """
    } else {
        """
	/* module load python2.7/Python-2.7.3_gnu446*/
	/* module load cutadapt/1.7.1_Python-2.7.3*/
     /restricted/projectnb/pulmseq/kkarri_netflow_test/trim_galore --paired --gzip $reads
        """
    }
}

process buildIndex
{
   publishDir "${params.outdir}/STAR", mode: 'copy'
    input:
    file params.genome
    file gtf
    output:
    file "STARgenome" into STARgenomeIndex
    """
    module load star/2.4.2a
    mkdir STARgenome
    STAR --runThreadN ${task.cpus} \
         --runMode genomeGenerate \
         --genomeDir STARgenome \
         --genomeFastaFiles ${params.genome} \
         --sjdbGTFfile ${gtf} \
         --sjdbOverhang 100 \
    """
}

process star {
	tag "$reads"
/*	executor = 'sge'
  * 	clusterOptions = "-P ${PROJECT} -l h_rt=96:00:00 -l mem_total=5G -pe omp 12"
*  cpus 10
  *  memory '80GB'
  *  time  { 5.h * task.attempt }
  *  errorStrategy { task.exitStatus == 143 ? 'retry' : 'terminate' }
  *  maxRetries 3
  *  maxErrors '-1'
*/
    publishDir "${params.outdir}/STAR", mode: 'copy'

    input:
    file index
    file gtf
    file (reads:'*') from trimmed_reads

    output:
    set file('*Log.final.out'), file ('*.bam') into aligned
    file '*.out' into star_logs
    file '*SJ.out.tab'

    script:
    """
	module load star/2.4.2a
	f='$reads';f=(\$f);f=\${f[0]};f=\${f%.gz};f=\${f%.fastq};f=\${f%.fq};f=\${f%_val_1};f=\${f%_trimmed};f=\${f%_1};f=\${f%_R1}
	prefix=\$f
	STAR --genomeDir $index \\
        --sjdbGTFfile $gtf \\
        --readFilesIn $reads  \\
        --runThreadN ${task.cpus} \\
        --twopassMode Basic \\
        --outWigType bedGraph \\
        --outSAMtype BAM SortedByCoordinate \\
        --readFilesCommand zcat \\
        --outFileNamePrefix \$prefix
 """
}

process stringtieFPKM {
    tag "$bam_stringtieFPKM"

 /*  clusterOptions = "-P ${PROJECT} -l h_rt=96:00:00 -l mem_total=5G -pe omp 12"

  *  memory { 4.GB * task.attempt }
   * time { 2.h * task.attempt }
   * errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
   * maxRetries 3
   * maxErrors '-1'
*/
    publishDir "${params.outdir}/stringtieFPKM", mode: 'copy'

    input:
    file bam_stringtieFPKM
    file gtf from gtf

    output:
    file '*_transcripts.gtf'
    file '*.gene_abund.txt'
    file '*.cov_refs.gtf'
    stdout into stringtie_log

    script:
    """
   /restricted/projectnb/pulmseq/kkarri_netflow_test/stringtie/stringtie $bam_stringtieFPKM \\
        -o ${bam_stringtieFPKM}_transcripts.gtf \\
        -v \\
        -G $gtf \\
        -A ${bam_stringtieFPKM}.gene_abund.txt \\
        -C ${bam_stringtieFPKM}.cov_refs.gtf \\
        -e \\
        -b ${bam_stringtieFPKM}_ballgown
 echo "File name: $bam_stringtieFPKM Stringtie version "\$(stringtie --version)
    """
}
def num_bams
bam_count.count().subscribe{ num_bams = it }

process multiqc {

    publishDir "${params.outdir}/MultiQC", mode: 'copy'


    input:
    file ('fastqc/*') from fastqc_results.toList()
    file ('trimgalore/*') from trimgalore_results.toList()
	file ('star/*') from star_logs.toList() 
    file ('stringtie/*') from stringtie_log.toList() 

    output:
    file '*multiqc_report.html'
    file '*multiqc_data'

    script:
    """
	module load python/2.7.11
	module load multiqc/0.8
    multiqc -f /restricted/projectnb/pulmseq/kkarri_netflow_test/results_new/fastqc/*.zip /restricted/projectnb/pulmseq/kkarri_netflow_test/results_new/trim_galore/ /restricted/projectnb/pulmseq/kkarri_netflow_test/results_new/stringtieFPKM/ /restricted/projectnb/pulmseq/kkarri_netflow_test/results_new/STAR/

    """
}

