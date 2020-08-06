exomiser_input_suffix = '\\.vcf\\.gz'
exomiser_analysisMode = 'PASS_ONLY'
exomiser_path = '/exomiser'
exomiser_exe = 'exomiser-cli-12.1.0.jar'
exomiser_extra_args = "--spring.config.location=${exomiser_path}/"
input_path = '/home/pontikos/d/jingyu/3041-3053'
exomiser_hpos = 'HP:0000556'

test_sample = 'B240054.vcf.gz'

Input_ch = Channel.fromPath("${input_path}/*.vcf.gz").map {it -> [it.name - ~/${exomiser_input_suffix}$/, it]}

process 'Exomiser' {
    tag "$sampleId"
    label "exomiser"
    container params.exomiser_docker
    //publishDir 's3://phenopolis-nextflow/exomiser/2020-08-04'
    publishDir "."
    // takes 30 minutes on 1 core / 10G
    cpus 1
    memory "10 G"

    input:
        tuple val(sampleId), path(input_file) from Input_ch
    output:
        path("${sampleId}.exomiser.tar.gz") into Exomiser_ch

    """
    set -e
    # make analysis.yml
    output_folder=${sampleId}-exomiser
    sed 's/samplePlaceholder/${input_file}/' ${exomiser_path}/exome-analysis-template.yml | sed 's/analysisModePlaceholder/${exomiser_analysisMode}/' | sed 's/hpoPlaceholder/${exomiser_hpos}/' | sed "s/outputPlaceholder/\${output_folder}/" > analysis.yml
    cat analysis.yml
    # exomiser
    mkdir \${output_folder}
    java -Xms4g -Xmx8g -jar ${exomiser_path}/${exomiser_exe} --analysis analysis.yml ${exomiser_extra_args}
    # zip up
    tar -zcvf ${sampleId}.exomiser.tar.gz \${output_folder}/
    
    """
}
