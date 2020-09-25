exomiser_path = '/exomiser'
exomiser_exe = 'exomiser-cli-12.1.0.jar'
exomiser_extra_args = "--spring.config.location=${exomiser_path}/"

Channel
  .fromPath(params.input_table)
  .splitCsv(header:['sampleId', 'read1', 'read2', 'hpo'], sep:',')
  .map{ row -> tuple(row.sampleId, row.read1, row.read2, row.hpo) }
  .set{ Input_ch }

process 'Exomiser' {
    tag "$sampleId"
    label "exomiser"
    container params.exomiser_docker
    // takes 30 minutes on 1 core / 10G
    cpus 1
    memory "10 G"

    input:
        tuple val(sampleId), val(read1), val(read2), val(hpo) from Input_ch
    output:
        path("exomiser.tar.gz") into Exomiser_ch

    script:
        exomiser_hpos = hpo.replaceAll(';',',')
    """
    set -e
    source s3.bash

    aws_profile="${params.s3_deposit_profile}"
    downloads=()
    for filename in ${params.input_filename} ${params.input_filename}.tbi; do
        cloudFile=${params.s3_deposit}/${sampleId}/\$filename
        downloads+=("nxf_s3_retry nxf_s3_download \$cloudFile \$filename")
    done
    nxf_parallel "\${downloads[@]}"

    # make analysis.yml
    output_folder=${sampleId}-exomiser
    sed 's/samplePlaceholder/${params.input_filename}/' ${exomiser_path}/exome-analysis-template.yml | sed 's/analysisModePlaceholder/${params.exomiser_analysisMode}/' | sed 's/hpoPlaceholder/${exomiser_hpos}/' | sed "s/outputPlaceholder/\${output_folder}/" > analysis.yml
    cat analysis.yml
    # exomiser
    mkdir \${output_folder}
    java -Xms4g -Xmx8g -jar ${exomiser_path}/${exomiser_exe} --analysis analysis.yml ${exomiser_extra_args}
    # zip up
    tar -zcvf  exomiser.tar.gz \${output_folder}/

    aws_profile="${params.s3_deposit_profile}"
    uploads=()
    uploads+=("nxf_s3_retry nxf_s3_upload exomiser.tar.gz ${params.s3_deposit}/${sampleId}")
    nxf_parallel "\${uploads[@]}"
    
    """
}
