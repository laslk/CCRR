#!/bin/bash

set -e
set -o pipefail


echo "********************"
echo "******Run SvABA*****"
echo "********************"


while getopts ":n:t:r:p:i:b:" opt; do
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
    i)
      id="$OPTARG"
      ;;
    p)
      threads="$OPTARG"
      echo "threads: $threads"
      ;;
    b)
      refbuild="$OPTARG"
      echo "build: $refbuild"
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

ref_dir=$(dirname "$refdata")
ref_base=$(basename "$refdata")


cd "$ref_dir"

if [ ! -f "${ref_base}.amb" ] || [ ! -f "${ref_base}.ann" ] || [ ! -f "${ref_base}.bwt" ] || [ ! -f "${ref_base}.pac" ] || [ ! -f "${ref_base}.sa" ]; then
  echo "BWA index files not found, generating..."
  bwa index "$ref_base"
fi

if [ ! -f "${ref_base}.fai" ]; then
  echo "FAI index file not found, generating..."
  samtools faidx "$ref_base"
fi


if [ "$refbuild" == "hg19" ]; then
  known_indels="${DATABASE_DIR}/database/dbsnp/hg19_v0_Homo_sapiens_assembly19.known_indels.vcf"
fi

if [ "$refbuild" == "hg38" ]; then
  known_indels="${DATABASE_DIR}/database/dbsnp/hg38_v0_Homo_sapiens_assembly38.known_indels.vcf"
fi

if [[ -d "${WORK_DIR}/${id}/svaba" ]]; then
  rm -rf ${WORK_DIR}/${id}/svaba
fi
mkdir ${WORK_DIR}/${id}/svaba
cd ${WORK_DIR}/${id}/svaba

svaba run -t $tumor -n $normal -p $threads -D $known_indels -a svaba -G $refdata || { rm -rf ${WORK_DIR}/${id}/svaba ; echo "Error during svaba"; exit 1; }


echo "********************"
echo "*****SvABA done*****"
echo "********************"
