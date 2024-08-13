#!/bin/bash

set -e  
set -o pipefail  

soreca="${TOOL_DIR}/soreca/bin/soreca"


echo "********************"
echo "*****Run Soreca*****"
echo "********************"


while getopts ":n:t:r:i:" opt; do
  case $opt in
    n)
      normal="$OPTARG"
	  echo "normal bam: $normal"
      ;;
    t)
      tumor="$OPTARG"
      echo "tumor bam: $tumor"
      ;;
    r)
      ref="$OPTARG"
      echo "reference: $ref"
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


if [[ -d "${WORK_DIR}/${id}/soreca" ]]; then
  rm -rf  ${WORK_DIR}/${id}/soreca
fi
mkdir ${WORK_DIR}/${id}/soreca
cd ${WORK_DIR}/${id}/soreca

$soreca enspan -i "$tumor" -o tumor -build "$ref" || { echo "Error during tumor enspan"; exit 1; }
$soreca enspan -i "$normal" -o normal -build "$ref" || { echo "Error during normal enspan"; exit 1; }
$soreca unsnarl -inT tumor_enspan.txt -inN normal_enspan.txt -o soreca -build "$ref" || { echo "Error during unsnarl"; exit 1; }

rm -f normal_enspan.bam      
rm -f normal_enspan.txt   
rm -f tumor_enspan.bam      
rm -f tumor_enspan.txt
rm -f normal_enspan.bam.bai  
rm -f tumor_enspan.bam.bai

echo "********************"
echo "****Soreca done*****"
echo "********************"

