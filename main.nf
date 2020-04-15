#!/usr/bin/env nextflow
/*
NXF_VER=20.01.0 nextflow run aws_align.nf -w s3://phenopolis-nextflow/nextflow-work -with-report report.html
refer to https://github.com/CRG-CNAG/CalliNGS-NF
and https://app.terra.bio/#workspaces/help-gatk/Germline-SNPs-Indels-GATK4-b37/workflows (wgs)
and https://app.terra.bio/#workspaces/help-gatk/Exome-Analysis-Pipeline (wes)
version >= 20.01.0

Exome/b37 only for now. But will adapt for genome/38
*/

/*
 * Define default parameters
 * Note that do_excessHet only applies to samples of unrelated samples with no consanguinity.
 * re: https://gatkforums.broadinstitute.org/gatk/discussion/23216/how-to-filter-variants-either-with-vqsr-or-by-hard-filtering
 * ExcessHet filtering applies only to callsets with a large number of samples, e.g. hundreds of unrelated samples. Small cohorts should not trigger ExcessHet filtering as values should remain small. Note cohorts of consanguinous samples will inflate ExcessHet, and it is possible to limit the annotation to founders for such cohorts by providing a pedigree file during variant calling.
 * ExcessHet for small number of samples is usually small and would still be safe to do this..
 * However for consanguinous or many related samples, ExcessHet inflates, and it might filter out good variants.
 */

/*
 *  PART 1: Alignment

 *  process 1A fastqc.
 Note that trimming is not needed for exome / genome. 
 ref1: https://gatkforums.broadinstitute.org/gatk/discussion/2957/read-trimming
 ref2: https://www.biostars.org/p/212136/
 */

 
params.gatk_docker= 'broadinstitute/gatk:4.1.4.1'
params.picard_docker= 'broadinstitute/picard:2.22.3'
params.bamOutput = 's3://phenopolis-test'
//params.gatk_bundle = 's3://broad-references/hg19/v0'
params.gatk_bundle = 's3://phenopolis-nextflow/gatk_bundle'
params.gatk_options = "-Djava.io.tmpdir=."
params.human_ref_base = 'Homo_sapiens_assembly19'
params.human_ref = "${params.human_ref_base}.fasta"
params.BQSR_known_sites = 'dbsnp_138.b37.vcf.gz Mills_and_1000G_gold_standard.indels.b37.vcf.gz Homo_sapiens_assembly19.known_indels_20120518.vcf'
params.input = 's3://phenopolis-nextflow/fastq/*_{1,2}.fastq.gz'
params.result = 's3://phenopolis-nextflow'
params.build = 'b37'

log.info """\
C A L L I N G S  -  N F    v 1.0
================================
genome   : $params.build
"""
Reads_ch = Channel.fromFilePairs(params.input)
//human_refs = Channel.fromPath("${params.gatk_bundle}/Homo_sapiens_assembly19.fasta*")
//reads = "s3://${params.cloudPath}${params.input}"
//GATK_bundle = Channel.value(file("${params.gatk_bundle}/*"))
Human_ref_ch = Channel.value(file("${params.gatk_bundle}/${params.human_ref}*"))

/*
 * process 1A align
 */
process '1A_align' {
  container 'job-definition://nf-biocontainers-bwa-v0-7-17-3-deb_cv1'
  containerOptions '-m 25g'
  input:
    tuple val(sampleId), path(reads) from Reads_ch
    file '*' from Human_ref_ch
  
  output:
    tuple sampleId, path("raw.bam") into Bwa_bam_ch
  
  """
  bwa mem -K 100000000 -v 3 -t 14 -Y ${params.human_ref} ${reads} 2> >(tee bwa.stderr.log >&2) \
    | \
		samtools view -1 - > raw.bam
  """
}
/*
 * process 1B mark duplicate
 */
process '1B_mark_duplicate' {
  memory '20 GB'
  container params.gatk_docker
  containerOptions '-m 14g'
  input:
    tuple val(sampleId), path(input_bam) from Bwa_bam_ch
  
  output:
    tuple sampleId, path("unsorted.duplicates_marked.bam") into Mark_duplicate_ch

  
  """
  gatk --java-options \" -Dsamjdk.compression_level=5 -Xms6000m\" \
    MarkDuplicates \
      --INPUT ${input_bam} \
      --OUTPUT unsorted.duplicates_marked.bam \
      --METRICS_FILE duplicate.metrics \
      --VALIDATION_STRINGENCY SILENT \
      --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
      --ASSUME_SORT_ORDER \"queryname\" \
      --CREATE_MD5_FILE true
  """
}

/*
 * process 1C add group
 */
process '1C_addGroup' {
  container params.picard_docker

  input:
    tuple val(sampleId), path(input_bam) from Mark_duplicate_ch
  output:
    tuple sampleId, path("group.bam") into Group_bam_ch

  """
  java -jar /usr/picard/picard.jar AddOrReplaceReadGroups \
    INPUT=${input_bam} \
    OUTPUT=group.bam \
    RGLB=unknown \
    RGPL=Illumina \
    RGPU=1 \
    RGSM=unknown \
    VALIDATION_STRINGENCY=SILENT
  """
}

/*
 * process 1D Sort BAM file by coordinate order, and fix readgroup
 */
process '1D_sort' {
  memory '30 G'
  container params.gatk_docker
  containerOptions '-m 20g'
  input:
    tuple val(sampleId), path(input_bam) from Group_bam_ch
  
  output:
    tuple sampleId, path("sorted.bam"), path("sorted.bam.bai") into Sorted_bam_ch
  
  """
  gatk --java-options \"${params.gatk_options} -Dsamjdk.compression_level=5 -Xms6000m\" \
    SortSam \
    --TMP_DIR=. \
    --INPUT ${input_bam} \
    --OUTPUT sorted.bam \
    --SORT_ORDER \"coordinate\" \
    --CREATE_INDEX true \
    --CREATE_MD5_FILE true
  samtools index sorted.bam
  """
}

// Duplicate Sorted_bam_ch
Sorted_bam_ch.into {Sorted_bam_ch1; Sorted_bam_ch2}
// Get gatk bundle files
GATK_ch = Channel.value(file("${params.gatk_bundle}/*"))
/*
 * process 1E BQSR
 */
process '1E_BQSR' {
  memory '12 G'
  container params.gatk_docker
  input:
    tuple val(sampleId), path(input_bam), path(input_bam_index) from Sorted_bam_ch1
    file '*' from GATK_ch
  output:
    path("recal_data.csv") into bqsr_report_ch
  
  """
  # -L can be given for parallelism
  # make known sites a string
  known_site_list=(${params.BQSR_known_sites})
  known_site_list=\$(printf " --known-sites %s" "\${known_site_list[@]}")
  known_site_list=\${known_site_list:1}
  gatk --java-options \"${params.gatk_options} -Xms4000m\" \
    BaseRecalibrator \
      -R ${params.human_ref} \
      -I ${input_bam} \
      --use-original-qualities \
      -O recal_data.csv \
      \$known_site_list 
  """
}

/*
 * process 1F Apply BQSR
 */
Human_ref_withdict_ch = Channel.value(file("${params.gatk_bundle}/${params.human_ref_base}*"))
process '1F_ApplyBQSR' {
  memory '10 G'
  container params.gatk_docker
  publishDir "${params.result}/bam", mode: 'copy'
  input:
    tuple val(sampleId), path(input_bam), path(input_bam_index) from Sorted_bam_ch2
    path(recal_file) from bqsr_report_ch
    file '*' from Human_ref_withdict_ch
  
  output:
    tuple sampleId, path("${sampleId}.bqsr.bam"), path("${sampleId}.bqsr.bam.bai") into BQSR_bam_ch
  
  """
  gatk --java-options \"${params.gatk_options} -Xms3000m\" \
    ApplyBQSR \
      -R ${params.human_ref} \
      -I ${input_bam} \
      -O ${sampleId}.bqsr.bam \
      -bqsr ${recal_file} \
      --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 \
      --add-output-sam-program-record \
      --create-output-bam-md5 \
      --use-original-qualities
  samtools index ${sampleId}.bqsr.bam
  """
}
