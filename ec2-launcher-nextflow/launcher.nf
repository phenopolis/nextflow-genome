to_dos = params.do.tokenize(',').join('|')
now = new Date().format( 'yyyy-MM-dd' )

process 'Stage' {
    cpus 1
    output:
        val true into Stage_ch
    """
    WORK_DIR=`mktemp -d -p "."`
    git -C \$WORK_DIR clone ${params.git_repo}
    ${params.aws_exe} s3 --profile nextflow cp \${WORK_DIR}/${params.git_repo_name} ${params.s3_repo} --recursive
    rm -rf \${WORK_DIR}
    """
}
process 'Align' {
    label 'aws'
    cpus 1
    memory '6 G'
    container params.nextflow_docker
    input:
        val stageDone from Stage_ch
    output:
        tuple val(true), path("report-align.html")  into Align_ch
        path("nextflow*.log") into AlignLog_ch

    """
    curl -s https://get.nextflow.io | bash
    if [[ algin =~ ^(${to_dos})\$ ]]; then
        /home/ec2-user/miniconda/bin/aws s3 cp ${params.s3_repo}/align . --recursive
        NXF_VER=${params.nextflow_version} ./nextflow -log nextflow-align.log run align.nf -bucket-dir ${params.bucket_dir}/${now}/align -with-report report-align.html --input_table ${params.input_table}
    else
        touch report-align.html
        touch nextlfow-align.log
    fi
    """
}


process 'variantCall' {
    label 'aws'
    cpus 1
    memory '6 G'
    container params.nextflow_docker
    echo true
    input:
        tuple val(alignDone), path("*") from Align_ch
        path("*") from AlignLog_ch
    output:
        tuple val(true), path("*.html") into VariantCall_ch
        path("nextflow*.log") into VariantCallLog_ch

    """
    curl -s https://get.nextflow.io | bash
    if [[ variantCall =~ ^(${to_dos})\$ ]]; then
        /home/ec2-user/miniconda/bin/aws s3 cp ${params.s3_repo}/variantCall . --recursive
        /home/ec2-user/miniconda/bin/aws s3 cp ${params.input_table} ./input_table.csv
        echo ${params.input_table}
        NXF_VER=${params.nextflow_version} ./nextflow -log nextflow-variantCall.log run variantCall.nf -bucket-dir ${params.bucket_dir}/${now}/variantCall -with-report report-variantCall.html --input_table ./input_table.csv
    else
        touch report-variantCall.html
        touch nextflow-variantCall.log
    fi
    """
}

process 'annotation' {
    label 'aws'
    cpus 1
    memory '6 G'
    container params.nextflow_docker
    publishDir '.'
    input:
        tuple val(variantCallDone), path("*") from VariantCall_ch
        path("*") from VariantCallLog_ch
    output:
        path("*.html") into Annotation_ch
        path("nextflow*.log") into AnnotationLog_ch

    """
    curl -s https://get.nextflow.io | bash
    if [[ annotation =~ ^(${to_dos})\$ ]]; then
        /home/ec2-user/miniconda/bin/aws s3 cp ${params.s3_repo}/annotation . --recursive
        NXF_VER=${params.nextflow_version} ./nextflow -log nextflow-annotation.log run annotation.nf -bucket-dir ${params.bucket_dir}/${now}/annotation -with-report report-annotation.html --input_table ${params.input_table}
    else
        touch report-annotation.html
        touch nextflow-annotation.log
    fi
    """
}
