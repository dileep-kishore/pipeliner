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
params.reads_file = params.genomes[params.genome].sample_reads_file
params.paired = params.genomes[params.genome].paired
params.outdir = params.genomes[params.genome].outdir
params.starindex= params.genome ? params.genomes[params.genome].star ?: false : false
params.bowtieindex= params.genome ? params.genomes[params.genome].bowtie ?: false : false
params.project = params.genome ? params.genomes[params.genome].project ?: false : false
params.aligner = params.genome ? params.genomes[params.genome].aligner ?: false : false
params.email = params.genome ? params.genomes[params.genome].email ?: false : false

// Read and Map reads with samples using csv file
def read_tuples = []
if (params.paired) {
    Channel
        .fromPath(params.reads_file)
        .splitCsv(header:true)
        .subscribe {row ->
            read_tuples.add(new Tuple(
                row.Sample_Name,
                new Tuple(
                    new File(row.Read1).absolutePath,
                    new File(row.Read2).absolutePath
                    )
                )
            )
        }
}
else {
    Channel
        .fromPath(params.reads_file)
        .splitCsv(header:true)
        .subscribe {row ->
            read_tuples.add(new Tuple(
                row.Sample_Name,
                new Tuple(
                    new File(row.Read).absolutePath
                    )
                )
            )
        }
}

println(read_tuples)

// Aligner options
if (params.aligner != 'star' && params.aligner != 'bowtie'){
    exit 1, "Invalid aligner option: ${params.aligner}. Valid options: 'star','bowtie'"
}



// Validate inputs
if (params.starindex && params.aligner == 'star' && false){
    starindex = Channel
        .fromPath(params.starindex)
        .ifEmpty { exit 1, "STAR index not found: ${params.starindex}" }
        .toList()
}
else if (params.bowtieindex && params.aligner == 'bowtie' && false){
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
if (params.gtf && false){
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
log.info "FASTA       : ${params.fasta}"
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
    .ifEmpty { error "File ${params.reads_file} not parsed properly" }
    .into { read_files_fastqc; read_files_trimming }

/* PREPROCESSING - Build STAR index */

if(params.aligner == 'star' && !params.starindex && fasta && false)
{
println("index does not exits")
process makeSTARindex {
        tag fasta
        publishDir "ref_genome/",  mode: 'copy'


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
Channel
    .fromFilePairs('reads/ggal_*{1,2}.fq')
    .into {temp_fasta_reads}

process fastqc {
    tag "$sampleid"
    publishDir "${params.outdir}/$sampleid/fastqc", mode: 'copy'

    input:
    // set val(name),file(reads) from read_files_fastqc
    set sampleid, reads from read_files_fastqc
    // set sampleid, file(reads) from temp_fasta_reads

    output:
    set sampleid, '*_fastqc.{zip,html}' into fastqc_results

    script:
    if (params.paired) {
        """
        module load fastqc/0.11.5
        fastqc -o . -q ${reads[0]} ${reads[1]}
        """
    }
    else {
        """
        /restricted/projectnb/pulmseq/kkarri_netflow_test/FastQC/fastqc -q ${reads[0]}
        """
    }
}




/*
 * STEP 2 - Trim Galore!
 */
process trim_galore {
 tag "$sampleid"

    publishDir "${params.outdir}/$sampleid/trim_galore", mode: 'copy'

    input:
    // set val(name),file(reads) from read_files_trimming
    set sampleid, reads from read_files_trimming

    output:
    set sampleid, '*fq.gz' into trimmed_reads
    set sampleid, '*trimming_report.txt' into trimgalore_results

    script:
    if(params.paired) {
        println("pairedend")
        """
        module load trim_galore/0.4.2
        trim_galore --paired --gzip ${reads[0]} ${reads[1]}
        """
    } else {
        println("singleend")
        """
        /* module load python2.7/Python-2.7.3_gnu446*/
        /* module load cutadapt/1.7.1_Python-2.7.3*/
        /restricted/projectnb/pulmseq/kkarri_netflow_test/trim_galore --gzip ${reads[0]}
        """
    }
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
