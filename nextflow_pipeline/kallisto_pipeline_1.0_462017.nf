#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
/*
 * SET UP CONFIGURATION VARIABLES
 */


// Pipeline version
version = 0.9

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
params.transcriptome = params.genome ? params.genomes[params.genome].transcriptome ?: false : false
params.kallistoindex = params.genome ? params.genomes[params.genome].kallistoindex ?: false : false
params.email = params.genome ? params.genomes[params.genome].email ?: false : false
params.refgenome = params.genome ? params.genomes[params.genome].refgenome ?: false : false
params.bed = params.genome ? params.genomes[ params.genome ].bed ?: false : false
params.saveReference = true

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
                    row.Read1,
                    row.Read2
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
                    row.Read
                    )
                )
            )
        }
}


// Aligner options

if (params.aligner != 'star' && params.aligner != 'bowtie' && params.aligner != 'kallisto'){
    exit 1, "Invalid aligner option: ${params.aligner}. Valid options: 'star','bowtie','kallisto'"
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

else if(params.kallistoindex && params.aligner == 'kallisto'){

   kallistoindex = Channel
        .fromPath(params.kallistoindex)
        .ifEmpty { exit 1, "Kallisto index not found: ${params.kallistoindex}" }
        .toList()
}


else if ( params.fasta ){
    fasta = file(params.fasta)
        if( !fasta.exists() ) exit 1, "Fasta file not found: ${params.fasta}"
}

else if ( ( params.aligner == 'bowtie' && !params.bowtieindex ) && !params.fasta ){
        exit 1, "No reference genome specified!"
}

if (params.transcriptome) {

        Channel
        .fromPath(params.transcriptome)
        .ifEmpty { exit 1, "Transcriptome file not found: ${params.transcriptome}" }
        .toList()
        .into { make_kallistoindex; transcriptome_kallisto  }
}


if (params.gtf){
    Channel
        .fromPath(params.gtf)
        .ifEmpty { exit 1, "GTF annotation file not found: ${params.gtf}" }
        .toList()
        .into { gtf_makeSTARindex; gtf_bowtieindex;gtf_star;gtf_bowtie;gtf_stringtieFPKM;  }
}

if( params.bed ){
    bed = Channel
        .fromPath(params.bed)
        .ifEmpty { exit 1, "BED annotation file not found: ${params.bed}" }
        .toList()

}
fasta = file(params.fasta)
gtf   = file(params.gtf)
transcriptome = file(params.transcriptome)
bed = file(params.bed)



log.info " PIPELINER RNAseq Pipeline v${version}"
log.info "===================================="
log.info "Project       : ${params.project}"
log.info "Reads         : ${params.reads}"
log.info "Genome        : ${params.genome}"
log.info "FASTA         : ${params.fasta}"
log.info "Transcriptome : ${params.transcriptome}"
log.info "Annotation    : ${params.gtf}"
log.info "Output dir    : ${params.outdir}"

        if(params.aligner == 'star')
        {
        log.info "Aligner        : STAR"
                if (params.starindex)
                        log.info "STAR Index: ${params.starindex}"

        }
        else if (params.aligner == 'bowtie')
        {
        log.info "Aligner        : Bowtie"
                if (params.bowtieindex)
                        log.info "Bowtie Index   : ${params.bowtieindex}"

        }

        else if (params.aligner == 'kallisto')
        {
        log.info "Pseudo Aligner        : Kallisto"
                if (params.kallistoindex)
                        log.info "Kallisto Index   : ${params.kallistoindex}"

        }

log.info "Current home  : $HOME"
log.info "Current user  : $USER"
log.info "Current path  : $PWD"
log.info "===================================="


Channel
    .from(read_tuples)
    .ifEmpty { error "File ${params.reads_file} not parsed properly" }
    .into { read_files_fastqc; read_files_trimming }


/* PREPROCESSING - Build STAR index */

if(params.aligner == 'star' && !params.starindex &&fasta)

{ 
process makeSTARindex {

		publishDir path: "${params.outdir}/reference_genome", saveAs: { params.saveReference ? it : null }, mode: 'copy'


		tag fasta
		cache 'deep'
                
		input:
                file fasta from fasta
                file gtf from gtf_makeSTARindex
		
                output:
                file "star" into starindex

                script:
		
                """
		module load star/2.4.2a
		mkdir star
                STAR \\
                --runMode genomeGenerate \\
                --runThreadN ${task.cpus} \\
                --sjdbGTFfile $gtf \\
                --sjdbOverhang 149 \\
                --genomeDir star/ \\
                --genomeFastaFiles $fasta

                """
                }

}

else if(params.aligner == 'kallisto' && !params.kallistoindex)


{               process Kallisto_index {

                tag transcriptome
                publishDir path: "${params.outdir}/reference_genome/kallisto", saveAs: { params.saveReference ? it : null }, mode: 'copy'


                input:
                file transcriptome_file from make_kallistoindex

                output:
                file "index" into index

                script:
                """
		module load  kallisto/
		
                kallisto index -i index ${transcriptome_file}

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
		cache 'deep'
	    	publishDir "${params.outdir}/$sampleid/fastqc", mode: 'copy'
		  executor = 'sge'
  		clusterOptions = "-P ${params.project} -l h_rt=96:00:00 -l mem_total=5G -pe omp 12"
	

   		input:
    		set sampleid, reads from read_files_fastqc

    		output:
    		set sampleid, '*_fastqc.{zip,html}' into fastqc_results

    		script:
    		if (params.paired) 
		{
        
			"""
        		/restricted/projectnb/pulmseq/kkarri_netflow_test/FastQC/fastqc -o . -q ${reads[0]} ${reads[1]}

        		"""
    		}
    	
		else 
		{
        
			"""
        		/restricted/projectnb/pulmseq/kkarri_netflow_test/FastQC/fastqc -o . -q ${reads[0]}
        		"""
    		}
	}


/*
 * STEP 2 - Trim Galore!
 */

process trim_galore {
 	tag "$sampleid"
	cache 'deep'
	  executor = 'sge'
  	clusterOptions = "-P ${params.project} -l h_rt=96:00:00 -l mem_total=5G -pe omp 12"

	
	publishDir "${params.outdir}/$sampleid/trim_galore", mode: 'copy'

    	input:
    	set sampleid, reads from read_files_trimming

    	output:
    	set sampleid, '*fq.gz' into trimmed_reads
    	set sampleid, '*trimming_report.txt' into trimgalore_results

    	script:
    	if(params.paired) 
		{
        
			println("pairedend")
        	
			"""
        		module load python2.7/Python-2.7.3_gnu446
        		module load cutadapt/1.7.1_Python-2.7.3
        		/restricted/projectnb/pulmseq/kkarri_netflow_test/trim_galore --paired --gzip ${reads[0]} ${reads[1]}
        		"""
    		} 
	else 
		{
        		println("singleend")
        		"""
        		/* module load python2.7/Python-2.7.3_gnu446*/
        		/* module load cutadapt/1.7.1_Python-2.7.3*/
        		/restricted/projectnb/pulmseq/kkarri_netflow_test/trim_galore --gzip ${reads[0]}
        		"""
    		}
	}


if (params.aligner == "kallisto")
{

process kallisto_mapping {



        publishDir "${params.outdir}/$sampleid/kallisto", mode: 'copy'

        input:

        file index from index.first()
        set sampleid, file (reads:'*') from trimmed_reads

        output:
        file "kallisto_${sampleid}" into kallisto_out_dirs

        script:

        if( params.paired ) {

                """
                module load anaconda2/4.3.0
                module load kallisto/0.43.0


                kallisto quant --paired -i index -o kallisto_${sampleid} ${reads}
                """
                        }
        else    {
                """
                module load anaconda2/4.3.0
                module load kallisto/0.43.0


                kallisto quant --single -i ${transcriptome_index} -o kallisto_${sampleid} ${reads}
                """
                        }

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
