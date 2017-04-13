#!/usr/bin/env nextflow
// vim: syntax=groovy

version = 0.2

params.paired = true
params.param_file = file('paired_params.csv')
test_file = file('test.sh');

if (params.paired) {
    Channel
        .fromPath(params.param_file)
        .splitCsv(header: true)
        .map {row -> [row.Sample_Name, [row.Read1, row.Read2]]}
        .into {foo; bar}
} else {
    Channel
        .fromPath(params.param_file)
        .splitCsv(header: true)
        .map {row -> [row.Sample_Name, [row.Read]]}
        .into {foo; bar}
}

process tester {
    echo true
    tag "$site"
    publishDir "$site", mode: "copy"
    input:
    set val(site), reads from foo
    output:
    set site, "*.txt" into other_reads
    file "*.txt" into oop_reads
    script:
    if (params.paired) {
        """
        echo ${reads[0]}
        head -n 20 ${reads[0]} ${reads[1]}> ${site}_paired_head.txt
        tail -n 20 ${reads[0]} ${reads[1]}> ${site}_paired_tail.txt
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

// oop_reads
//     .flatten()
//     .map {f -> f.path}
//     .collectFile(name:'dragon1.txt', newLine:true)
//     .into {sample_results}

process makeafolder {
    echo true
    publishDir "final", mode: "copy"
    input:
    // set site, file (ht_file:'*') from oop_reads
    // val temp_list from other_reads.collect {it[1]}
    val temp_list from oop_reads.collect()
    output:
    file "*.txt" into temp_reads
    script:
    String temp_string = temp_list .flatten() .join(' ')
    // String temp_list = temp
    //                     .findAll {e -> e.class == ArrayList}
    //                     .flatten()
    //                     .join(' ')
    """
    echo $temp_string
    cat $temp_string > dragon2.txt
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
