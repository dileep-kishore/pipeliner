#!/usr/bin/env nextflow
// vim: syntax=groovy

version = 0.1

params.param_file = "params.txt"
def samples = []
def file_names = []
Channel
    .fromPath(params.param_file)
    .splitCsv(header:true)
    // .flatten()
    .subscribe {row -> file_names.add(new File(row.FastQ_files).absolutePath); samples.add(row.sample_name)}

def sample_files = []
for(i=0;i<file_names.size;i++) {sample_files.add([samples[i],file_names[i]])}
def sample_map = sample_files.groupBy {it[0]}
                             .collectEntries { key, value -> new Tuple( key, value*.getAt(1) ) }
def file_pairs = []
for(i in sample_map) {
    file_pairs.add(new Tuple(i.key, i.value[0], i.value[1]))
    // file_pairs.add(new Tuple(i.key, i.value))
}
// sample_map.each{ key, values -> file_pairs.add(new Tuple(key, values)) }

// file_pairs.each { elem -> println(elem[1].absolutePath) }

Channel
    .from(file_pairs)
    // .subscribe{ println(it) }
    .into {foo; bar}

params.paired = true
if (params.paired) {
    params.read_files = "ggal_*_{1,2}.fq"
    Channel
        .fromFilePairs(params.read_files)
        // .subscribe { println(it)}
        .set{rna_reads}
    }
else {
    // params.read_files = "ggal_*_1.fq"
    params.read_files = "ggal_*_{1,2}.fq"
    Channel
        .fromPath(params.read_files)
        .set{rna_reads}
    }

rna_reads.into {single_reads; paired_reads}

if (!params.paired)
{
    process output_single {
        tag "$reads"
        echo true
        input:
        file reads from single_reads
        output:
        file reads into op_reads
        script:
        """
        echo $reads
        """
    }
}


if (params.paired)
{
    process output_paired {
        echo true
        tag "$site"
        input:
        set site, file(reads) from paired_reads
        output:
        file reads into op_reads
        script:
        """
        """
    }
}

process tester {
    echo true
    tag "$site"
    publishDir "$site", mode: "copy"
    input:
    set site, read1, read2 from bar
    output:
    file "*.txt" into oop_reads
    script:
    """
    echo $site
    echo $read1 $read2
    head -n 20 $read1 > temp.txt
    """
}

// op_reads.into {test1, test2}

// process output_check {
//     echo true
//     input:
//     file op_reads
//     script:
//     """
//     echo $op_reads
//     """
// }
