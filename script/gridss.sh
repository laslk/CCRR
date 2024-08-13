#!/bin/bash

set -e
set -o pipefail


gridss_database="${DATABASE_DIR}/database/gridss_database/"


echo "********************"
echo "*****Run gridss*****"
echo "********************"


while getopts ":n:t:r:p:s:b:j:1:2:i:" opt; do
  case $opt in
    n)
      normal=$(realpath "$OPTARG")
      echo "normal bam: $normal"
      ;;
    t)
      tumor=$(realpath "$OPTARG")
      echo "tumor bam: $tumor"
      ;;
    r)
      refdata=$(realpath "$OPTARG")
      echo "reference: $refdata"
      ;;
    p)
      threads="$OPTARG"
      echo "threads: $threads"
      ;;
    i)
      id="$OPTARG"
      ;;
    s)
      skip="$OPTARG"
      echo "If using BWA, skipping the soft clip realignment step:$skip"
      ;;
    b)
      ref="$OPTARG"
      echo "reference build: $ref"
      ;;
    j)
      jvmheap="$OPTARG"
      echo "jvm heap:$jvmheap G"
      ;;
    1)
      name_n="$OPTARG"
      echo "normal name:$name_n"
      ;;
    2)
      name_t="$OPTARG"
      echo "tumor name:$name_t"
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done



if [[ -d "${WORK_DIR}/${id}/gridss" ]]; then
  rm -rf "${WORK_DIR}/${id}/gridss"
fi

mkdir ${WORK_DIR}/${id}/gridss
cd ${WORK_DIR}/${id}/gridss

if [ "$ref" == "hg19" ];then
  if [ "$skip" == "True" ]; then
    gridss -t $threads -r $refdata -o gridss.vcf --skipsoftcliprealignment -b ${TOOL_DIR}/gridss/example/ENCFF001TDO.bed $normal $tumor --jvmheap ${jvmheap}g || { rm -rf ${WORK_DIR}/${id}/gridss;echo "Error during gridss"; exit 1; }
  elif [ "$skip" == "False" ];then
    gridss -t $threads -r $refdata -o gridss.vcf -b ${TOOL_DIR}/gridss/example/ENCFF001TDO.bed $normal $tumor --jvmheap ${jvmheap}g || { rm -rf ${WORK_DIR}/${id}/gridss ; echo "Error during gridss"; exit 1; }
  fi
elif [ "$ref" == "hg38" ];then
  if [ "$skip" == "True" ]; then
    gridss -t $threads -r $refdata -o gridss.vcf --skipsoftcliprealignment -b ${TOOL_DIR}/gridss/example/ENCFF356LFX.bed $normal $tumor --jvmheap ${jvmheap}g || { rm -rf ${WORK_DIR}/${id}/gridss;echo "Error during gridss"; exit 1; }
  elif [ "$skip" == "False" ];then
    gridss -t $threads -r $refdata -o gridss.vcf -b ${TOOL_DIR}/gridss/example/ENCFF356LFX.bed $normal $tumor --jvmheap ${jvmheap}g || { rm -rf ${WORK_DIR}/${id}/gridss ; echo "Error during gridss"; exit 1; }
  fi
fi

rm -rf *.working
rm -rf *.bam

if [[ ! -d "${WORK_DIR}/${id}/gridss/gripss_output" ]]; then
  mkdir ${WORK_DIR}/${id}/gridss/gripss_output
fi

if [ "$ref" == "hg19" ];then
  java -Xmx${jvmheap}G -jar /opt/conda/envs/main/share/hmftools-gripss-2.3.2-0/gripss.jar -sample $name_t -reference $name_n -ref_genome_version 37 -ref_genome $refdata -pon_sgl_file ${gridss_database}37/sv/sgl_pon.37.bed.gz -pon_sv_file ${gridss_database}37/sv/sv_pon.37.bedpe.gz -known_hotspot_file ${gridss_database}37/sv/known_fusions.37.bedpe -repeat_mask_file ${gridss_database}37/sv/repeat_mask_data.37.fa.gz -vcf gridss.vcf -output_dir gripss_output  || { echo "Error during gripss"; exit 1; }
elif [ "$ref" == "hg38" ];then
  java -Xmx${jvmheap}G -jar /opt/conda/envs/main/share/hmftools-gripss-2.3.2-0/gripss.jar -sample $name_t -reference $name_n -ref_genome_version 38 -ref_genome $refdata -pon_sgl_file ${gridss_database}38/sv/sgl_pon.38.bed.gz -pon_sv_file ${gridss_database}38/sv/sv_pon.38.bedpe.gz -known_hotspot_file ${gridss_database}38/sv/known_fusions.38.bedpe -repeat_mask_file ${gridss_database}38/sv/repeat_mask_data.38.fa.gz -vcf gridss.vcf -output_dir gripss_output  || { echo "Error during gripss"; exit 1; }
fi


echo "********************"
echo "*****gridss done****"
echo "********************"