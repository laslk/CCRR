#!/bin/bash

echo "******************"
echo "***CTLPScanner****"
echo "******************"

while getopts ":b:1:c:o:m:" opt; do
  case $opt in
    b)
      refbuild="$OPTARG"
      ;;
    1)
      tumor_id="$OPTARG"
      ;;
    c)
      cn="$OPTARG"
      ;;
    o)
      wd="$OPTARG"
      ;;
    m)
      mt="$OPTARG"
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

if [[ -d $wd ]]; then
  rm -rf $wd
fi

mkdir $wd
cd $wd

python ${SCRIPT_DIR}/format_conv/cn2ctlpscanner.py $cn ${wd}/ctlp_input.txt --sample-name $tumor_id || { echo "Error during CTLPscanner"; exit 1; }
echo "${wd}/ctlp_input.txt" > ${wd}/inputlist || { echo "Error during CTLPscanner"; exit 1; }
Rscript ${TOOL_DIR}/ctlpscanner/CTLP_detection.R $refbuild $wd wgs
mv ${wd}/ctlpscanner_results.txt ${wd}/CTLPRegion.txt
cp ${TOOL_DIR}/ctlpscanner/drawSVG.pl ${wd}/drawSVG.pl
perl ${wd}/drawSVG.pl

rm -rf ${wd}/drawSVG.pl
rm -rf inputlist
rm -rf ctlp_input.txt

