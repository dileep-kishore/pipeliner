#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Pipeline version
version = 0.1

// Configurable variables
params.genome = "/restricted/projectnb/pulmseq/kkarri_netflow_test/genome/genome.index.fa"
params.gtf  = "/restricted/projectnb/pulmseq/kkarri_netflow_test/genome/genome.bed.gff"
params.reads = "/restricted/projectnb/pulmseq/kkarri_netflow_test/Data/ggal/*{_1,_2}.fq"
params.outdir = "/restricted/projectnb/pulmseq/kkarri_netflow_test/results_new"
params.index = "/restricted/projectnb/pulmseq/kkarri_netflow_test/STARgenome"
/*PROJECT = params.project*/



log.info "===================================="
log.info " RNAseq Pipeline v${version}"
log.info "===================================="
log.info "Reads       : ${params.reads}"
log.info "Genome      : ${params.genome}"
log.info "Annotation   : ${params.gtf}"
log.info "Output dir   : ${params.outdir}"
log.info "Index       : ${params.index}"
log.info "===================================="

// Validate inputs
index = file(params.index)
gtf   = file(params.gtf)
if( !index.exists() ) exit 1, "Missing STAR index: $index"
if( !gtf.exists() )   exit 2, "Missing GTF annotation: $gtf"

/*
 * Create a channel for input read files
 */
Channel
    .fromPath( params.reads )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .set { read_files }

read_files.into {read_files_fastqc; read_files_trimming}


/*
 * STEP 1 - FastQC
 */
process fastqc {

/*executor = 'sge'
*	clusterOptions = "-P ${PROJECT} -l h_rt=24:00:00 -l mem_total=5G"
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


/*
 * STEP 2 - Trim Galore!
 */
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

/* Step 3  STAR Build Index */

process buildIndex {

  
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


/* Step 4  STAR Alignment */


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

// Function that checks the alignment rate of the STAR output
// and returns true if the alignment passed and otherwise false
def check_log(logs) {
    def percent_aligned = 0;
    logs.eachLine { line ->
        if ((matcher = line =~ /Uniquely mapped reads %\s*\|\s*([\d\.]+)%/)) {
            percent_aligned = matcher[0][1]
        }
    }
    if(percent_aligned.toFloat() <='10'.toFloat() ){
        println "#################### VERY POOR ALIGNMENT RATE ONLY ${percent_aligned}%! FOR ${logs}"
        false
    } else {
        println "Passed aligment with ${percent_aligned}%! FOR ${logs}"
        true
    }
}

// Filter removes all 'aligned' channels that fail the check
aligned
    .filter { logs, bams -> check_log(logs) }
    .flatMap {  logs, bams -> bams }
    .set { SPLIT_BAMS }
SPLIT_BAMS.into { bam_count; bam_rseqc; bam_preseq; bam_markduplicates; bam_featurecounts; bam_stringtieFPKM }

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



	
workflow.onComplete {
        println ( workflow.success ? "CONGRATULATIONS !!!!! Your pipeline executed successfully :) !!" : "Oops .. something went wrong" )
}
