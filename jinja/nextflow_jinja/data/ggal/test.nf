#!/usr/bin/env nextflow
// vim: syntax=groovy

version = 0.2

params.paired = false
params.param_file = "single_params.csv"
test_file = file("test.sh");

def file_tuple = []
if (params.paired) {
    Channel
        .fromPath(params.param_file)
        .splitCsv(header:true)
        .subscribe {row ->
            file_tuple.add(new Tuple(
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
        .fromPath(params.param_file)
        .splitCsv(header:true)
        .subscribe {row ->
            file_tuple.add(new Tuple(
                row.Sample_Name,
                new Tuple(
                    new File(row.Read).absolutePath
                    )
                )
            )
        }
    }


Channel
    .from(file_tuple)
    // .subscribe{ println(it) }
    // .flatMap { a,b,c -> [a, c]}
    // .collate( 2 )
    // .subscribe { println(it) }
    .into {foo; bar}

process tester {
    echo true
    tag "$site"
    publishDir "$site", mode: "copy"
    input:
    set site, reads from foo
    output:
    set site, "*.txt" into oop_reads
    script:
    if (params.paired) {
        """
        echo ${reads[0]}
        head -n 20 ${reads[0]} ${reads[1]}> paired_head.txt
        tail -n 20 ${reads[0]} ${reads[1]}> paired_tail.txt
        """
    }
    else {
        """
        echo ${reads[0]}
        head -n 20 ${reads[0]} > single_head.txt
        tail -n 20 ${reads[0]} > single_tail.txt
        """
    }
}

process makeafolder {
    echo true
    tag "$site"
    publishDir "$site", mode: "copy"
    input:
    set site, file (ht_file:'*') from oop_reads
    output:
    set site, "*.txt" into temp_reads
    script:
    def temp_list = [1,2,3]
    """
    echo $ht_file > temp.txt
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
