to_dos = params.do.tokenize(',').join('|')
now = new Date().format( 'yyyy-MM-dd' )

process 'Stage' {
    cpus 1
    echo true
    output:
        val true into Stage_ch
    """
    WORK_DIR=`mktemp -d -p "/tmp"`
    git -C \$WORK_DIR clone ${params.git_repo}
    ${params.aws_exe} s3 cp \${WORK_DIR}/${params.git_repo_name} ${params.s3_repo} --recursive
    rm -rf \${WORK_DIR}
    """
}
process 'Align' {
    label 'aws'
    cpus 1
    memory '6 G'
    container params.nextflow_docker
    echo true
    input:
        val stageDone from Stage_ch
    output:
        tuple val(true), path("report-align.html")  into Align_ch

    """
    if [[ algin =~ ^(${to_dos})\$ ]]; then
        /home/ec2-user/miniconda/bin/aws s3 cp ${params.s3_repo}/align . --recursive
        ${params.nextflow_cliPath} run align.nf -bucket-dir ${params.bucket_dir}/${now}/align -with-report report-align.html
    else
        touch report-align.html
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
        val alignDone from Align_ch
    output:
        tuple val(true), path("report-variantCall.html") into VariantCall_ch

    """
    if [[ variantCall =~ ^(${to_dos})\$ ]]; then
        /home/ec2-user/miniconda/bin/aws s3 cp ${params.s3_repo}/variantCall . --recursive
        ${params.nextflow_cliPath} run align.nf -bucket-dir ${params.bucket_dir}/${now}/variantCall -with-report report-variantCall.html
    else
        touch report-variantCall.html
    fi
    """
}

process 'annotation' {
    label 'aws'
    cpus 1
    memory '6 G'
    container params.nextflow_docker
    echo true
    input:
        val variantCallDone from VariantCall_ch
    output:
        tuple val(true), path("report-annotation.html") into Annotation_ch

    """
    if [[ annotation =~ ^(${to_dos})\$ ]]; then
        /home/ec2-user/miniconda/bin/aws s3 cp ${params.s3_repo}/annotation . --recursive
        ${params.nextflow_cliPath} run align.nf -bucket-dir ${params.bucket_dir}/${now}/annotation -with-report report-annotation.html
    else
        touch report-annotation.html
    fi
    """
}
