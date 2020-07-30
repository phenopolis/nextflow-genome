padd_bed() {
  input_bed=$1
  output_bed=$2
  padding=$3
  [ -f ${output_bed} ] && rm ${output_bed}
  while read line; do
    chrom=$(echo $line | awk '{print $1}')
    start=$(echo $line | awk '{print $2}')
    start=$(( ${start}-$padding ))
    start=$(( ${start} > 0 ? ${start} : 0 ))
    end=$(echo $line | awk '{print $3}')
    end=$(( ${end}+$padding ))
    echo -e "${chrom}\t${start}\t${end}" >> ${output_bed}
  done<${input_bed}
}

#padd_bed log/exome_target_00.bed log/output.bed 100
if_has_record() {
  infile=$1
  echo $(zcat ${infile} | grep -m 1 -v "#" -c)
}
