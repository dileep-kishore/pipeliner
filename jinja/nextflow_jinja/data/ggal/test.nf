#!/usr/bin/env nextflow
// vim: syntax=groovy

version = 0.1

params.param_file = "params.txt"
def file_pairs = []
Channel
    .fromPath(params.param_file)
    .splitCsv(header:true)
    .subscribe {row ->
        file_pairs.add(new Tuple(
            row.Sample_Name,
            new File(row.Read1).absolutePath,
            new File(row.Read2).absolutePath
            )
        )
    }

Channel
    .from(file_pairs)
    // .subscribe{ println(it) }
    // .flatMap { a,b,c -> [a, c]}
    // .collate( 2 )
    // .subscribe { println(it) }
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
    // publishDir "$site", mode: "copy"
    input:
    // set site, file(reads:'*') from bar
    set site, read1, read2 from foo
    output:
    set site, "*.txt" into oop_reads
    script:
    """
    head -n 20 $read1 > head.txt
    tail -n 20 $read2 > tail.txt
    """
}

process makeafolder {
    echo true
    tag "$site"
    publishDir "$site", mode: "copy"
    input:
    set site, file (reads:'*') from oop_reads
    output:
    set site, "*.txt" into temp_reads
    script:
    """
    echo $reads > temp.txt
    echo $site
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
