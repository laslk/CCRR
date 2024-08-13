#!/bin/bash

set -e
set -o pipefail


echo "********************"
echo "*****run cnvkit*****"
echo "********************"


while getopts ":n:t:r:i:b:u:l:p:" opt; do
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
      ;;
    u)
      cellularity="$OPTARG"
      echo "purity: $cellularity"
      ;;
    l)
      ploidy="$OPTARG"
      echo "ploidy: $ploidy"
      ;;
    b)
      ref="$OPTARG"
      echo "reference build: $ref"
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


if [[ -d "${WORK_DIR}/${id}/cnvkit" ]]; then
  rm -rf ${WORK_DIR}/${id}/cnvkit
fi
mkdir ${WORK_DIR}/${id}/cnvkit
cd ${WORK_DIR}/${id}/cnvkit

cnvkit.py batch $tumor -n $normal -m wgs -p $threads -f $refdata --annotate ${DATABASE_DIR}/database/refFlat/${ref}/refFlat.csv || { rm -rf ${WORK_DIR}/${id}/cnvkit;echo "Error during cnvkit"; exit 1; } 

cns_file=$(ls *.cns | grep -E '^[^.]+\.cns$')
if [ -z "$cns_file" ]; then
    echo "Error: No .cns file found."
    exit 1
else
    if [[ "$cellularity" != "null" && "$ploidy" != "null" ]]; then
        cnvkit.py call $cns_file --ploidy $(printf "%.0f" "$ploidy") --purity $cellularity -o recalibrated.call.cns || { rm -rf ${WORK_DIR}/${id}/cnvkit;echo "Error during cnvkit recalibrate"; exit 1; } 
    else
        cnvkit.py call $cns_file -o recalibrated.call.cns || { rm -rf ${WORK_DIR}/${id}/cnvkit;echo "Error during cnvkit recalibrate"; exit 1; } 
    fi
    python ${SCRIPT_DIR}/format_conv/cnvkit2bed.py recalibrated.call.cns ${id}.cnvkit.bed
fi

echo "*****************"
echo "***cnvkit done***"
echo "*****************"