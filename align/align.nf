#!/usr/bin/env nextflow
/*
NXF_VER=20.01.0 nextflow run align.nf -w s3://phenopolis-nextflow/nextflow-work --input_table <tab-delimited table> -with-report report.html
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


log.info """\
C A L L I N G S  -  N F    v 1.0
================================
genome   : $params.build
"""
// determines build
HG38 = ['hg38', 'GRCh38', 'grch38']
HG19 = ['hg19', 'GRCh37', 'b37', 'grch37']
human_ref_base = null
human_ref_path = null
gatk_bundle = null
if (HG38.contains(params.build)) {
  human_ref_base = params.human_ref_base_hg38
  human_ref_path = params.human_ref_path_hg38
  gatk_bundle = params.gatk_bundle_hg38
  BQSR_known_sites = params.BQSR_known_sites_hg38
} else if (HG19.contains(params.build)) {
  human_ref_path = params.human_ref_path_hg19
  human_ref_base = params.human_ref_base_hg19
  gatk_bundle = params.gatk_bundle_hg19
  BQSR_known_sites = params.BQSR_known_sites_hg19
} else {
  exit 1, "cannot understand build. Please provide either --build hg19, or --build hg38 (default is hg19)"
}
human_ref = "${human_ref_base}.fasta"
bwa_human_ref_bundle = "'${human_ref}' '${human_ref}.fai' '${human_ref}.pac' '${human_ref_base}.dict' '${human_ref}.amb' '${human_ref}.ann' '${human_ref}.bwt' '${human_ref}.sa'"
human_ref_bundle = "'${human_ref}' '${human_ref}.fai' '${human_ref}.pac' '${human_ref_base}.dict' '${human_ref}.amb' '${human_ref}.ann'"

// get input channel from input table
Channel
  .fromPath(params.input_table)
  .splitCsv(header:['sampleId', 'read1', 'read2'], sep:',')
  .map{ row -> tuple(row.sampleId, row.read1, row.read2) }
  .set{ Reads_ch }
/*
 * process 1A align
 */
process '1A_align' {
  // needs a way to throw an error if there's something wrong with bwa
  cpus 20
  memory '30 G'
  tag "$sampleId"
  label 'small_batch'
  container params.align_docker
  input:
    tuple val(sampleId), val(read1), val(read2) from Reads_ch
  
  output:
    tuple sampleId, path("raw.bam") into Bwa_bam_ch
  
  """
  set -e
  source s3.bash
  # copy human_ref
  aws_profile="${params.aws_ref_profile}"
  downloads=()
  for downfile in ${bwa_human_ref_bundle}
  do
    downloads+=("nxf_s3_retry nxf_s3_download ${human_ref_path}/\$downfile ./\$downfile")
  done
  # copy fastq
  aws_profile="${params.fastq_path_profile}"
  local_read1=${read1}
  local_read2=${read2}
  local_read1=\${local_read1##*/}
  local_read2=\${local_read2##*/}
  downloads+=("nxf_s3_retry nxf_s3_download ${read1} ./\${local_read1}")
  downloads+=("nxf_s3_retry nxf_s3_download ${read2} ./\${local_read2}")
  nxf_parallel "\${downloads[@]}"
  ls -alh

  bwa mem -K 100000000 -v 3 -t 20 -Y ${human_ref} \${local_read1} \${local_read2} 2> >(tee bwa.stderr.log >&2) \
    | \
  samtools view -1 - > raw.bam
  """
}

/*
 * process 1B mark duplicate
 */
process '1B_mark_duplicate' {
  tag "$sampleId"
  memory '20 GB'
  container params.gatk_docker
  label 'SSD'
  cpus 4
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
  tag "$sampleId"
  cpus 1
  memory '8 GB'
  container params.picard_docker
  label 'small_batch'

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
    RGSM=$sampleId \
    VALIDATION_STRINGENCY=SILENT
  """
}

/*
 * process 1D Sort BAM file by coordinate order, and fix readgroup
 */
process '1D_sort' {
  tag "$sampleId"
  memory '30 G'
  container params.gatk_docker
  label 'SSD'
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

/*
 * process 1E BQSR
 */
process '1E_BQSR' {
  tag "$sampleId"
  label 'SSD'
  memory '12 G'
  container params.gatk_docker
  input:
    tuple val(sampleId), path(input_bam), path(input_bam_index) from Sorted_bam_ch1
  output:
    path("recal_data.csv") into bqsr_report_ch
  
  """
  source s3.bash
  # copy known_sites
  aws_profile=""
  downloads=()
  for downfile in ${BQSR_known_sites}
  do
    downloads+=("nxf_s3_retry nxf_s3_download ${gatk_bundle}/\$downfile ./\$downfile")
    if [[ "\$downfile" == *gz ]]
    then
      downloads+=("nxf_s3_retry nxf_s3_download ${gatk_bundle}/\${downfile}.tbi ./\${downfile}.tbi")
    elif [[ "\$downfile" == *vcf ]]
    then
      downloads+=("nxf_s3_retry nxf_s3_download ${gatk_bundle}/\${downfile}.idx ./\${downfile}.idx")
    fi
  done
  nxf_parallel "\${downloads[@]}"
  # copy human_ref
  downloads=()
  aws_profile="${params.aws_ref_profile}"
  for downfile in ${human_ref_bundle}
  do
    downloads+=("nxf_s3_retry nxf_s3_download ${human_ref_path}/\$downfile ./\$downfile")
  done
  nxf_parallel "\${downloads[@]}"

  # -L can be given for parallelism
  # make known sites a string
  known_site_list=(${BQSR_known_sites})
  known_site_list=\$(printf " --known-sites %s" "\${known_site_list[@]}")
  known_site_list=\${known_site_list:1}
  gatk --java-options \"${params.gatk_options} -Xms4000m\" \
    BaseRecalibrator \
      -R ${human_ref} \
      -I ${input_bam} \
      --use-original-qualities \
      -O recal_data.csv \
      \$known_site_list 
  """
}

/*
 * process 1F Apply BQSR
 */

process '1F_ApplyBQSR' {
  tag "$sampleId"
  memory '10 G'
  label 'SSD'
  container params.gatk_docker
  //publishDir './bam'
  input:
    tuple val(sampleId), path(input_bam), path(input_bam_index) from Sorted_bam_ch2
    path(recal_file) from bqsr_report_ch
  
  output:
    tuple sampleId, path("bqsr.bam"), path("bqsr.bam.bai") into BQSR_bam_ch
  
  """
  source s3.bash
  # copy human_ref
  aws_profile="${params.aws_ref_profile}"
  downloads=()
  for downfile in ${human_ref_bundle}
  do
    downloads+=("nxf_s3_retry nxf_s3_download ${human_ref_path}/\$downfile ./\$downfile")
  done
  nxf_parallel "\${downloads[@]}"
  gatk --java-options \"${params.gatk_options} -Xms3000m\" \
    ApplyBQSR \
      -R ${human_ref} \
      -I ${input_bam} \
      -O bqsr.bam \
      -bqsr ${recal_file} \
      --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 \
      --add-output-sam-program-record \
      --create-output-bam-md5 \
      --use-original-qualities
  samtools index bqsr.bam

  aws_profile="${params.s3_deposit_profile}"
  uploads=()
  uploads+=("nxf_s3_retry nxf_s3_upload bqsr.bam ${params.s3_deposit}/${sampleId}")
  uploads+=("nxf_s3_retry nxf_s3_upload bqsr.bam.bai ${params.s3_deposit}/${sampleId}")
  nxf_parallel "\${uploads[@]}"
  """
}

