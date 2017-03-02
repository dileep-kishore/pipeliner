#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-

/*
 * SET UP CONFIGURATION VARIABLES
 */



// Pipeline version
version = 0.6

// Configurable variables
params.genome = 'ggal'
params.fasta = params.genome ? params.genomes[params.genome].fasta ?: false : false
params.gtf  = params.genome ? params.genomes[params.genome].gtf ?: false : false
params.reads = params.genomes[params.genome].reads
params.sample_files = params.genomes[params.genome].sample_reads_file
params.outdir = params.genomes[params.genome].outdir
params.starindex= params.genome ? params.genomes[params.genome].star ?: false : false
params.bowtieindex= params.genome ? params.genomes[params.genome].bowtie ?: false : false
params.project = params.genome ? params.genomes[params.genome].project ?: false : false
params.aligner = params.genome ? params.genomes[params.genome].aligner ?: false : false
params.email = params.genome ? params.genomes[params.genome].email ?: false : false
println(params.reads)

// Read and Map reads with samples using csv file
def read_tuples = []
Channel
    .fromPath(params.sample_files)
    .splitCsv(header:true)
    .subscribe {row ->
        read_tuples.add(new Tuple(
            row.Sample_Name,
            new File(row.Read1).absolutePath,
            new File(row.Read2).absolutePath
            )
        )
    }

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
    if (params.starindex)    log.info "STAR Index: ${params.starindex}"
} else if (params.aligner == 'bowtie') {
    log.info "Aligner        : Bowtie"
    if (params.bowtieindex)        log.info "Bowtie Index   : ${params.bowtieindex}"

}

log.info "Current home   : $HOME"
log.info "Current user   : $USER"
log.info "Current path   : $PWD"
log.info "===================================="


 /* Create a channel for input read files*/

// Channel
//     .fromFilePairs( params.reads, size:-1 )
//     .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
//      .set{read }
//     read.into {read_files_fastqc; read_files_trimming}

Channel
    .from(read_tuples)
    .ifEmpty { error "File ${params.sample_files} not parsed properly" }
    .into { read_files_fastqc; read_files_trimming }

/* PREPROCESSING - Build STAR index */

if(params.aligner == 'star' && !params.starindex && fasta)
{
println("index does not exits")
process makeSTARindex {
        tag fasta
        publishDir "~/restricted/projectnb/pulmseq/kkarri_netflow_test/reference_genome/",  mode: 'copy'


        input:
        file fasta from fasta
        file gtf from gtf_makeSTARindex

        output:
        file "Stargenome" into starindex


        script:
        """

        module load star/2.4.2a

        mkdir Stargenome
        STAR \\
            --runMode genomeGenerate \\
            --runThreadN ${task.cpus} \\
            --sjdbGTFfile $gtf \\
            --sjdbOverhang 149 \\
            --genomeDir /restricted/projectnb/pulmseq/kkarri_netflow_test/reference_genome/ \\
            --genomeFastaFiles $fasta
        """
    }

}

else { println( "index exits")}

if(params.aligner == 'bowtie' && !params.bowtieindex)

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



logParams(params, "nextflow_paramters.txt")

def logParams(p, n) {
  File file = new File(n)
  file.write "Parameter:\tValue\n"

  for(s in p) {
     file << "${s.key}:\t${s.value}\n"
  }
}







/*
 * STEP 1 - FastQC
 */
process fastqc {
     tag "$sampleid"
     executor = 'sge'
     clusterOptions = "-P ${params.project} -l h_rt=96:00:00 -l mem_total=5G -pe omp 12"


    publishDir "${params.outdir}/$sampleid/fastqc", mode: 'copy'

    input:
    // set val(name),file(reads) from read_files_fastqc
    set sampleid, read1, read2 from read_files_fastqc
    // file(reads:'*') from read_files_fastqc

    output:
    set sampleid, '*_fastqc.{zip,html}' into fastqc_results

    script:
    """
    /restricted/projectnb/pulmseq/kkarri_netflow_test/FastQC/fastqc -q $read1 $read2
    """
}




/*
 * STEP 2 - Trim Galore!
 */
process trim_galore {
 tag "$sampleid"
 executor = 'sge'
 clusterOptions = "-P ${params.project} -l h_rt=96:00:00 -l mem_total=5G -pe omp 12"

    publishDir "${params.outdir}/$sampleid/trim_galore", mode: 'copy'

    input:
	// set val(name),file(reads) from read_files_trimming
    set sampleid, read1, read2 from read_files_trimming
    //  file(reads:'*') from read_files_trimming

    output:
    set sampleid, '*fq.gz' into trimmed_reads
    set sampleid, '*trimming_report.txt' into trimgalore_results

    script:
    single = reads instanceof Path

    if(!single) {
                println("pairedend")
          """
        module load python2.7/Python-2.7.3_gnu446
        module load cutadapt/1.7.1_Python-2.7.3

     /restricted/projectnb/pulmseq/kkarri_netflow_test/trim_galore --paired --gzip $read1 $read2
        """
    } else {
        println("singleend")
        """
        /* module load python2.7/Python-2.7.3_gnu446*/
        /* module load cutadapt/1.7.1_Python-2.7.3*/

     /restricted/projectnb/pulmseq/kkarri_netflow_test/trim_galore --gzip $read1 $read2
        """
    }
}


if(params.aligner == 'star'){

process star {
        tag "$sampleid"
        executor = 'sge'
        clusterOptions = "-P ${params.project} -l h_rt=96:00:00 -l mem_total=5G -pe omp 12"




    publishDir "${params.outdir}/$sampleid/STAR", mode: 'copy'

    input:
    file index from starindex.first()
    file gtf from gtf_star.first()
    set sampleid, file (reads:'*') from trimmed_reads

    output:
    set sampleid, file('*Log.final.out'), file ('*.bam') into aligned
    set sampleid, '*.out' into alignment_logs
    set sampleid, '*SJ.out.tab' into alignment_tab

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



// TODO: Find out whether the bam files need to be aggregated
aligned
    .filter { sampleid, logs, bams -> check_log(logs) }
    .flatMap {  sampleid, logs, bams -> sampleid, bams }
    .collate( 2 )
    .set { SPLIT_BAMS }
SPLIT_BAMS.into { bam_count; bam_rseqc_bamstats; bam_rseqc_genecoverage; bam_rseqc_junc_annot; bam_preseq;  bam_stringtieFPKM }




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
    tag "$sampleid"

        executor = 'sge'
        clusterOptions = "-P ${params.project} -l h_rt=96:00:00 -l mem_total=5G -pe omp 12"



    publishDir "${params.outdir}/$sampleid/stringtieFPKM", mode: 'copy'

    input:
    set sampleid, bamfiles from bam_stringtieFPKM
    file gtf from gtf

    output:
    file '*_transcripts.gtf'
    file '*.gene_abund.txt'
    file '*.cov_refs.gtf'
    stdout into stringtie_log

    script:
    """
   /restricted/projectnb/pulmseq/kkarri_netflow_test/stringtie/stringtie $bamfiles \\
        -o ${bam_stringtieFPKM}_transcripts.gtf \\
        -v \\
        -G $gtf \\
        -A ${bam_stringtieFPKM}.gene_abund.txt \\
        -C ${bam_stringtieFPKM}.cov_refs.gtf \\
        -e \\
        -b ${bam_stringtieFPKM}_ballgown
 echo "File name: $bamfiles Stringtie version "\$(stringtie --version)
    """
}
def num_bams
bam_count.count().subscribe{ num_bams = it }


// BamStats for the bam file
process bam_stats{
    tag "$sampleid"

	executor = 'sge'
        clusterOptions = "-P ${params.project} -l h_rt=96:00:00 -l mem_total=5G -pe omp 12"

    publishDir "${params.outdir}/$sampleid/bam_stats", mode: 'copy'

    input:
    set sampleid, bamfiles from bam_rseqc_bamstats

    output:
    file '*info.txt' into bam_stats_results

    script:
    """
    module load python
    module load rseqc/2.6.4
    bam_stat.py -i $bamfiles > bam_stats_info.txt
    """
}

// GeneBodyCoverage for the reference (house-keeping) genes
process gene_body_coverage{
    tag "$sampleid"

	executor = 'sge'
        clusterOptions = "-P ${params.project} -l h_rt=96:00:00 -l mem_total=5G -pe omp 12"

    publishDir "${params.outdir}/$sampleid/gene_coverage", mode: 'copy'

    input:
    set sampleid, bamfiles from bam_rseqc_genecoverage
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
    geneBody_coverage.py -r $gtf -i $bamfiles -o gene_coverage
    """
}

process junction_annotation{
    tag "$sampleid"

	executor = 'sge'
        clusterOptions = "-P ${params.project} -l h_rt=96:00:00 -l mem_total=5G -pe omp 12"

    publishDir "${params.outdir}/$sampleid/junction_annotation", mode: 'copy'

    input:
    set sampleid, bamfiles from bam_rseqc_junc_annot
    file gtf from ref_gene_model

    output:
    file 'junc_annot*' into junction_annotation_results

    // gtf must be the reference gene model
    // index file for bam must be available

    script:
    """
    module load python
    module load rseqc/2.6.4
    junction_annotation.py -i $bamfiles -o junc_annot -r $gtf
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

    // TODO: Multi-sample MultiQC (will it work if -f ${params.outdir})
    script:
    """
         module load python/2.7.11
        module load multiqc/0.8

   multiqc -f ${params.outdir}
	"""
}


workflow.onComplete {
	println ( workflow.success ? "CONGRATULATIONS !!!!! Your pipeline executed successfully :) !!" : "Oops .. something went wrong" )
	def subject = 'My pipeline execution'
	def recipient = (params.email)

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
