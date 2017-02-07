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
params.starindex= params.genome ? params.genomes[params.genome].star ?: false : false
params.bowtieindex= params.genome ? params.genomes[params.genome].bowtie ?: false : false
params.project = params.genome ? params.genomes[params.genome].project ?: false : false
params.aligner = params.genome ? params.genomes[params.genome].aligner ?: false : false
params.email = params.genome ? params.genomes[params.genome].email ?: false : false



// Aligner options

if (params.aligner != 'star' && params.aligner != 'bowtie'){
    exit 1, "Invalid aligner option: ${params.aligner}. Valid options: 'star','bowtie'"
}



// Validate inputs

if (params.starindex && params.aligner == 'star'){
    starindex = Channel
        .fromPath(params.starindex)
        .ifEmpty { exit 1, "STAR index not found: ${params.starindex}" }
        .toList()
}

else if (params.bowtieindex && params.aligner == 'bowtie' ){
    bowtieindex = Channel
        .fromPath(params.bowtieindex)
        .ifEmpty { exit 1, "Bowtie index not found: ${params.bowtieindex}" }
        .toList()
}


else if ( params.fasta ){
    fasta = file(params.fasta)
    if( !fasta.exists() ) exit 1, "Fasta file not found: ${params.fasta}"
}
else if ( ( params.aligner == 'bowtie' && !params.bowtieindex ) && !params.fasta ){
    exit 1, "No reference genome specified!"
}

fasta = file(params.fasta)

gtf   = file(params.gtf)

if (params.gtf){
    Channel
        .fromPath(params.gtf)
        .ifEmpty { exit 1, "GTF annotation file not found: ${params.gtf}" }
        .toList()
        .into { gtf_makeSTARindex; gtf_bowtieindex;gtf_star;gtf_stringtieFPKM;  }
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
if(params.aligner == 'star'){
    log.info "Aligner        : STAR"
    if (params.star_index)    log.info "STAR Index     : ${params.star_index}"
} else if (params.aligner == 'bowtie') {
    log.info "Aligner        : Bowtie"
    if (params.bowtieindex)        log.info "Bowtie Index   : ${params.bowtieindex}"
 
}

log.info "Current home   : $HOME"
log.info "Current user   : $USER"
log.info "Current path   : $PWD"
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


if(params.aligner == 'bowtie' && !params.bowtieindex && fasta)

{
process buildIndex {

   tag fasta
   publishDir "${params.outdir}/bowtie", mode: 'copy'

    input:
    file fasta from fasta

    output:
    file 'genome.index*' into bowtieindex
    file genome_file

    """
    module load bowtie2/2.2.9

    bowtie2-build -c ${genome_file} genome.index
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
SPLIT_BAMS.into { bam_count; bam_rseqc_bamstats; bam_rseqc_genecoverage; bam_rseqc_junc_annot; bam_preseq; bam_markduplicates; bam_featurecounts; bam_stringtieFPKM }


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


if(params.aligner == 'bowtie'){


process mapping {
    
	executor = 'sge'
 	clusterOptions = "-P ${params.project} -l h_rt=96:00:00 -l mem_total=5G -pe omp 12"

	tag "$bowtieindex"
    	publishDir "${params.outdir}/tophat_out", mode: 'copy'

   input:
	file fasta from fasta
	file bowtieindex from bowtie_index.first()
        file gtf from gtf_bowtieindex.first()
        file (reads:'*') from trimmed_reads
   

    output:
    set pair_id, "accepted_hits.bam" into bam
    set file('*Log.final.out'), file ('*.bam') into aligned


    """
        module load boost/1.58.0
        module load bowtie2/2.2.9
        module load tophat/2.1.1


    tophat2 -G $gtf $bowtieindex ${reads}
    """
}

}

process bam_stats{
    tag "$bam_rseqc_bamstats"

	executor = 'sge'
        clusterOptions = "-P ${params.project} -l h_rt=96:00:00 -l mem_total=5G -pe omp 12"

    publishDir "${params.outdir}/bam_stats", mode: 'copy'

    input:
    file bam_rseqc_bamstats

    output:
    file '*.txt' into bam_stats_results

    script:
    """
    module load python
    module load rseqc/2.6.4
    bam_stat.py -i $bam_rseqc_bamstats > bam_stats_info.txt
    """
}

process gene_body_coverage{
    tag "$bam_rseqc_genecoverage"

	executor = 'sge'
        clusterOptions = "-P ${params.project} -l h_rt=96:00:00 -l mem_total=5G -pe omp 12"

    publishDir "${params.outdir}/gene_coverage", mode: 'copy'

    input:
    file bam_rseqc_genecoverage
    file gtf from gtf

    output:
    file 'gene_coverage*' into gene_coverage_results
    stdout into gene_coverage_log

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
    file gtf from gtf

    output:
    file 'junc_annot*' into junction_annotation_results
    stdout into gene_coverage_log

    script:
    """
    module load python
    module load rseqc/2.6.4
    junction_annotation.py -i $bam_rseqc_junc_annot -o junc_annot -r $gtf
    """
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

process Kallisto_index {
    
	publishDir "${params.outdir}/kallisto", mode: 'copy'

        executor = 'sge'
        clusterOptions = "-P ${params.project} -l h_rt=96:00:00 -l mem_total=5G -pe omp 12"


    input:
    file Gallus_gallus.Galgal4.cdna.all.fa.gz  

    output:
    file transcripts.idx

    script:
    """
    #conda install -c bioconda kallisto=0.43.0 : to install in an environment
    
    kallisto index -i transcripts.idx Gallus_gallus.Galgal4.cdna.all.fa.gz
    """
}

process Kallisto_Quantification {
	executor = 'sge'
        clusterOptions = "-P ${params.project} -l h_rt=96:00:00 -l mem_total=5G -pe omp 12" 


    input:
    file transcripts.idx from kallisto
    file (reads:'*') from trimmed_reads
    
    output:
    file abundance.h5 into kallisto
    file abundance.tsv into kallisto
    file run_info.json into kallisto

    script:
    """
    kallisto quant -i {$params.outdir}/kallisto/ -o output {$reads} 
    """
}

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



	
workflow.onComplete {
        println ( workflow.success ? "CONGRATULATIONS !!!!! Your pipeline executed successfully :) !!" : "Oops .. something went wrong" )
    	def subject = 'My pipeline execution'
    	def recipient = {params.email}

    ['mail', '-s', subject, recipient].execute() << """

    Pipeline execution summary
    ---------------------------
    Completed at: ${workflow.complete}
    Duration    : ${workflow.duration}
    Success     : ${workflow.success}
    workDir     : ${workflow.workDir}
    exit status : ${workflow.exitStatus}
    Error report: ${workflow.errorReport ?: '-'}
    """
}

