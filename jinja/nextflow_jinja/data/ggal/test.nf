#!/usr/bin/env nextflow
// vim: syntax=groovy

version = 0.1

params.paired = false
if (params.paired) {
    params.read_files = "ggal_*_{1,2}.fq"
    Channel
        .fromFilePairs(params.read_files)
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
        echo $site
        echo $reads
        """
    }
}

process output_check {
    echo true
    input:
    file op_reads
    script:
    """
    echo $op_reads
    """
}
