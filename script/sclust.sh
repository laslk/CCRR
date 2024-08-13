#!/bin/bash

set -e
set -o pipefail

sclust="${TOOL_DIR}/sclust/bin/Sclust"
varscan2sclust="${SCRIPT_DIR}/format_conv/varscan2sclust.py"



echo "********************"
echo "*****Run Sclust*****"
echo "********************"


while getopts ":n:e:a:t:r:i:" opt; do
  case $opt in
    n)
      normal=$(realpath "$OPTARG")
      echo "normal bam: $normal"
      ;;
    e)
      test="$OPTARG"
      echo "test:$test"
      ;;
    i)
      id="$OPTARG"
      ;;
    a)
      alpha="$OPTARG"
      echo "alpha: $alpha"
      ;;
    t)
      tumor=$(realpath "$OPTARG")
      echo "tumor bam: $tumor"
      ;;
    r)
      ref="$OPTARG"
      echo "reference: $ref"
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


varscan_input="${WORK_DIR}/${id}/varscan/varscan.somatic.vcf"
if [[ -d "${WORK_DIR}/${id}/sclust" ]]; then
  rm -rf  ${WORK_DIR}/${id}/sclust
fi
mkdir ${WORK_DIR}/${id}/sclust
cd ${WORK_DIR}/${id}/sclust

for ((i=1; i<=24; i++))
do
    if ((i <= 22)); then
        chromosome="chr${i}"
    elif ((i == 23)); then
        chromosome="chrX"
    elif ((i == 24)); then
        chromosome="chrY"
    fi
    read_count=$(samtools idxstats $normal | grep "^${chromosome}\s" | awk '{print $3}')
    if [ "$read_count" -ne "0" ]; then
        $sclust bamprocess -t "$tumor" -n "$normal" -o sclust -part 2 -build $ref -r $chromosome || { echo "Error during sclust bamprocess"; exit 1; }
    fi
done
$sclust bamprocess -i sclust -o sclust || { echo "Error during sclust bamprocess generate a read-count file"; exit 1; }

python $varscan2sclust -i $varscan_input -o sclust_input_from_varscan.vcf

if [ "$test" != "test" ];then
  $sclust cn -rc sclust_rcount.txt -snp sclust_snps.txt -vcf sclust_input_from_varscan.vcf -alpha $alpha -o sclust_final -ns 1000 || { echo "Error during sclust cn"; exit 1; }
else
  cp ${DATABASE_DIR}/data/test/sclust_final_iCN.seg .
fi



echo "********************"
echo "****Sclust done*****"
echo "********************"
