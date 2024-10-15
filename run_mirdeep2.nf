#!/usr/bin/env nextflow


params.mature_miRNA_list = "$baseDir/Mature.fa"
params.reference_genome = "$baseDir/pb_b10_ill1.fasta"
params.fastq_files = "$baseDir/sample_info.csv"
params.multiqc = "$baseDir/multiqc"
params.outdir = "results"
params.fastqc_outdir = "results/fastqc_output"
params.mapper_config = "$baseDir/config.txt"
params.hairpin = "$baseDir/Hairpin.fa"
params.mirna_second = "$baseDir/Mature_ath.fa"



log.info """\
    M I R D E E P 2  P I P E L I N E
    ===================================
    mature_miRNA_list: ${params.mature_miRNA_list}
    reference_genome: ${params.reference_genome}
    fastq_files: ${params.fastq_files}
    outdir: ${params.outdir}
    ShortStack threads: ${params.threads}
    """
    .stripIndent(true)




process READ_FASTQS {
    debug true
    input:
    tuple val(sampleId), file(read1)

    script:
    """
    echo head --sample $sampleId --reads $read1 
    """
}



process FASTQC {
    tag "FastQC on $sample_id"
    publishDir params.fastqc_outdir, mode: 'copy'


    input:
     tuple val(sample_id), path(reads)
        

    output:
    path "fastqc_${sample_id}"

    script:
    """
    echo ${reads}
    echo ${sample_id}
    mkdir fastqc_${sample_id}
    fastqc -o fastqc_${sample_id} -f fastq -q ${reads}
    """


}


process MULTIQC {

    publishDir params.outdir, mode: 'copy'

    input:
    path '*'

    output:
    path 'multiqc_report.html'


    script:
    """
    multiqc .
    """

}












process DECOMPRESS_GZIP_FILES {
    tag "Gzip -d on $sample_id"
    publishDir params.outdir, mode: 'copy'


    input:
     tuple val(sample_id), path(reads)
        

    output:
    //path "*.fastq"
    tuple val(sample_id), path("*.fastq")

    script:
    """


    echo ${reads}
    echo ${sample_id}
    
    gzip -f -d ${reads} 
    """


}





process BUILD_BOWTIE_INDEX {
  tag "Build Bowtie index"
    //publishDir params.fastqc_outdir, mode: 'copy'


    input:
     path(reference_genome)
        

    output:
    path "*.ebwt"
    //tuple val(sample_id), path("*.fastq")

    script:
    """
    bowtie-build ${reference_genome} bowtie_index     
    """



}



process MAPPER{
     tag "Mapper.pl on $sample_id"
    publishDir params.outdir, mode: 'copy'


    input:
    // tuple val(sample_id), path(reads)
     path(config_file)
     path(bowtie_index)
    
    output:
    path "reads_collapsed.fa" , emit: reads_collapsed
    path "reads_collapsed_vs_genome.arf" , emit: arf

    script:
    """
    
    #bwt_idx=\$( basename bowtie_index.1.ebwt .1.ebwt)
    #mapper.pl ${config_file} -d -e -h -j   -l 18 -m -p \${bwt_idx} -s reads_collapsed.fa -t reads_collapsed_vs_genome.arf -v -q -o 12
    mapper.pl ${config_file} -d -e -h -j   -l 18 -m -p bowtie_index -s reads_collapsed.fa -t reads_collapsed_vs_genome.arf -v -q -o 2
    """
}


process MIRDEEP2{
     tag "MIRDEEP2 on $sample_id"
    publishDir params.outdir, mode: 'copy'


    input:
     path(reads_collapsed)
     path(arf)
     path(genome)
     path(mature)
     path(mature_second)
     path(hairpin)

    output:
    path "report.log" 

    script:
    """

    miRDeep2.pl ${reads_collapsed} ${genome} ${arf} ${mature} ${mature_second} ${hairpin} -g -1 -c 2> report.log
    """
}


workflow {


    samples_ch = Channel.fromPath(params.fastq_files)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample_id, file(row.read1)) }

    
    collected_samples_ch = Channel.fromPath(params.fastq_files)
        .splitCsv(header: true)
        .map { row -> tuple(file(row.read1)) }
        .collect()
        .view()



   // fastqc_ch = FASTQC(samples_ch)
   // MULTIQC(fastqc_ch.collect())
   
   // MIRTRACE(collected_samples_ch)
  
   //fastq_ch =  DECOMPRESS_GZIP_FILES(samples_ch)

    bowtie_index_ch = BUILD_BOWTIE_INDEX(params.reference_genome)
    mapper_ch = MAPPER(params.mapper_config,bowtie_index_ch )
    mirdeep2_ch =  MIRDEEP2(MAPPER.out.reads_collapsed, MAPPER.out.arf, params.reference_genome, params.mature_miRNA_list, params.mirna_second, params.hairpin )
}


workflow.onComplete {
    log.info ( workflow.success ? "\nDone" : "Oops .. something went wrong" )

}