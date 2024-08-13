#!/bin/bash

set -e
set -o pipefail

echo "********************"
echo "*Run AmpliconSuite**"
echo "********************"


while getopts ":n:t:c:r:d:p:i:" opt; do
  case $opt in
    n)
      normal="$OPTARG"
      echo "normal bam: $normal":
      ;;
    t)
      tumor="$OPTARG"
      echo "tumor bam: $tumor"
      ;;
    c)
      cnv="$OPTARG"
      echo "cnv: $cnv"
      ;;
    r)
      ref="$OPTARG"
      echo "reference: $ref"
      ;;
    d)
      threads="$OPTARG"
      echo "threads: $threads"
      ;;
    p)
      prefix="$OPTARG"
      echo "prefix: $prefix"
      ;;
    i)
      id="$OPTARG"
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

if [[ -z "$normal" || -z "$tumor" || -z "$ref" || -z "$cnv" || -z "$prefix" || -z "$threads" ]]; then
  echo "One or more required arguments are missing." >&2
  exit 1
fi

if [[ ! -d "${WORK_DIR}/${id}/complex/AA" ]]; then
  mkdir ${WORK_DIR}/${id}/complex/AA
fi
cd ${WORK_DIR}/${id}/complex/AA

python ${SCRIPT_DIR}/format_conv/cn2AAinput.py -input $cnv -output ${WORK_DIR}/${id}/complex/AA/cnvinput.bed || { echo "Error during AA"; exit 1; }
AmpliconSuite-pipeline.py -s $prefix -t $threads --cnv_bed ${WORK_DIR}/${id}/complex/AA/cnvinput.bed \
    --bam $tumor --ref $ref --normal_bam $normal --run_AA --run_AC || { echo "Error during AA"; exit 1; }


echo "********************"
echo "******AA done*******"
echo "********************"
