#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-

/*
 * SET UP CONFIGURATION VARIABLES
 */



// Pipeline version
version = 0.2

// Configurable variables
params.genome = 'ggal'
params.fasta = params.genome ? params.genomes[params.genome].fasta ?: false : false
params.gtf  = params.genome ? params.genomes[params.genome].gtf ?: false : false
params.reads = params.genomes[params.genome].reads
params.outdir = params.genomes[params.genome].outdir
params.index= params.genome ? params.genomes[params.genome].star ?: false : false
params.project = params.genome ? params.genomes[params.genome].project ?: false : false




// Choose aligner
params.aligner = 'star'
if (params.aligner != 'star'){
    exit 1, "Invalid aligner option: ${params.aligner}. Valid options: 'star'"
}



// Validate inputs
if( params.index && params.aligner == 'star'){
    index = Channel
        .fromPath(params.index)
        .ifEmpty { exit 1, "STAR index not found: ${params.index}" }
        .toList()
}

if ( params.fasta ){
    fasta = file(params.fasta)
    if( !fasta.exists() ) exit 1, "Fasta file not found: ${params.fasta}"
}


gtf   = file(params.gtf)

if( params.gtf ){
    Channel
        .fromPath(params.gtf)
        .ifEmpty { exit 1, "GTF annotation file not found: ${params.gtf}" }
        .toList()
        .into { gtf_makeSTARindex; gtf_star;gtf_stringtieFPKM }
}



log.info "===================================="
log.info " PIPELINER RNAseq Pipeline v${version}"
log.info "===================================="
log.info "Project       : ${params.project}"
log.info "Reads       : ${params.reads}"
log.info "Genome      : ${params.genome}"
log.info "FASTA	      : ${params.fasta}"
log.info "Annotation   : ${params.gtf}"
log.info "Output dir   : ${params.outdir}"
log.info "Index       : ${params.index}"
log.info "===================================="


 /* Create a channel for input read files*/

Channel
    .fromPath( params.reads )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .set { read_files }

read_files.into { read_files_fastqc; read_files_trimming}


/* PREPROCESSING - Build STAR index */

if(params.aligner == 'star' && !params.star_index && fasta){
    process makeSTARindex {
        tag fasta
        publishDir "~/restricted/projectnb/pulmseq/kkarri_netflow_test/", mode: 'copy'

        input:
        file fasta from fasta
        file gtf from gtf_makeSTARindex

        output:
        file "STARgenome" into star_index
        
        script:
        """
        module load star/2.4.2a

        mkdir STARgenome
        STAR \\
            --runMode genomeGenerate \\
            --runThreadN ${task.cpus} \\
            --sjdbGTFfile $gtf \\
            --sjdbOverhang 149 \\
            --genomeDir STARgenome/ \\
            --genomeFastaFiles $fasta
        """
    }
}



/*
 * STEP 1 - FastQC
 */
process fastqc {

     executor = 'sge'
     clusterOptions = "-P ${params.project} -l h_rt=96:00:00 -l mem_total=5G -pe omp 12"


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


 executor = 'sge'
 clusterOptions = "-P ${params.project} -l h_rt=96:00:00 -l mem_total=5G -pe omp 12"

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



/* Step 3  STAR Alignment */

if(params.aligner == 'star'){

process star {
	tag "$reads"
	executor = 'sge'
        clusterOptions = "-P ${params.project} -l h_rt=96:00:00 -l mem_total=5G -pe omp 12"




    publishDir "${params.outdir}/STAR", mode: 'copy'

    input:
    file index from star_index.first()
    file gtf from gtf_star.first()
    file (reads:'*') from trimmed_reads

    output:
    set file('*Log.final.out'), file ('*.bam') into aligned
    file '*.out' into alignment_logs
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

aligned
    .filter { logs, bams -> check_log(logs) }
    .flatMap {  logs, bams -> bams }
    .set { SPLIT_BAMS }
SPLIT_BAMS.into { bam_count; bam_rseqc; bam_preseq; bam_markduplicates; bam_featurecounts; bam_stringtieFPKM }


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







process stringtieFPKM {
    tag "$bam_stringtieFPKM"

	executor = 'sge'
        clusterOptions = "-P ${params.project} -l h_rt=96:00:00 -l mem_total=5G -pe omp 12"



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



/*
 * STEP 11 MultiQC
 */
process multiqc {
     
	executor = 'sge'
        clusterOptions = "-P ${params.project} -l h_rt=96:00:00 -l mem_total=5G -pe omp 12"

	publishDir "${params.outdir}/MultiQC", mode: 'copy'

    input:
    file ('fastqc/*') from fastqc_results.flatten().toList()
    file ('trimgalore/*') from trimgalore_results.flatten().toList()
    file ('alignment/*') from alignment_logs.flatten().toList()
    file ('stringtie/*') from stringtie_log.flatten().toList()
    

    output:
    file "*multiqc_report.html"
    file "*multiqc_data"

    script:
    """
	 module load python/2.7.11
        module load multiqc/0.8
	
   multiqc -f ${params.outdir}/fastqc/ ${params.outdir}/STAR/ ${params.outdir}/trim_galore/ ${params.outdir}/stringtieFPKM/
    """
}



	
workflow.onComplete {
        println ( workflow.success ? "CONGRATULATIONS !!!!! Your pipeline executed successfully :) !!" : "Oops .. something went wrong" )
}




