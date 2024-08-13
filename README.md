### CCRR: Complex Chromosomal Rearrangements Resolver
Complex Chromosomal Rearrangements Resolver (CCRR), a user-friendly workflow for analyzing complex chromosomal rearrangements in tumors. CCRR provides a simple, versatile, flexible, and reproducible solution for detecting somatic complex rearrangement events in tumors. Based on Docker containerization technology, CCRR allows users to freely choose and install tools for detecting SV and CNV. It automatically installs the tools, configures dependencies, and uses a merge algorithm to obtain accurate consensus SVs and CNVs, enabling the execution of various complex rearrangement callers. With a single command, users can obtain a comprehensive overview of complex rearrangements in the tumor genome
##### 1. Installation
###### Download
```
wget -O ccrr1.1.zip http://life-bioinfo.tpddns.cn:13362/download_file/ccrr1.1.zip
unzip -q ccrr1.1.zip
```
###### Prepare for installation using `install.py`
```
python install.py -sequenza -manta -delly -svaba -gridss \
                    -lumpy -soreca -purple -sclust -cnvkit \
                    -ref 'hg19'
```
This script will automatically download the dependency data for the tools you selected and build the Dockerfile.
```
  -sequenza             use sequenza for cellularity and ploidy

  -manta                use manta for sv
  -delly                use delly for sv and cn
  -svaba                use svaba for sv
  -gridss               use gridss for sv
  -lumpy                use lumpy for sv
  -soreca               use soreca for sv
  -purple               use purple for cn
  -sclust               use sclust for cn
  -cnvkit               use cnvkit for cn

  -ref                  hg19 or 'hg19&hg38'
```
You should select at least one SV tool and one CN tool.If you want to run in the fastest way, you can use only Delly to obtain SV and CN.

###### Obtain licenses for Mosek and Gurobi 
Gurobi: Apply for a WLS Compute Server license and store it in the same directory as the Dockerfile, named `gurobi.lic`. For more information, visit www.gurobi.com
Mosek: Obtain a Mosek license and store it in the same directory as the Dockerfile, named `mosek.lic`.  For more information, visit www.mosek.com
###### Build Docker image 
```
docker build --pull --rm --build-arg GITHUB_PAT=[GITHUB_PAT] \
        --build-arg SCRIPT_DIR="/home/0.script" --build-arg TOOL_DIR="/home/1.tools" \
        --build-arg DATABASE_DIR="/home/2.share" --build-arg WORK_DIR="/home/3.wd" \
        -f Dockerfile -t ccrr:v1 .
```
To ensure the installation proceeds correctly, you need to provide a GitHub token `GITHUB_PAT`.
You can specify these four parameters: `SCRIPT_DIR` for the script directory, default is `/home/0.script`; `TOOL_DIR` for the tool directory, default is `/home/1.tools`; `DATABASE_DIR` for the database and mounted shared directory, default is `/home/2.share`; `WORK_DIR` for the work directory, default is `/home/3.wd`.
###### Run a container 
```
docker run -v $(pwd)/share:/home/2.share -v $(pwd)/wd:/home/3.wd -d -it --name ccrr ccrr:v1
docker exec -it ccrr /bin/bash
```
This command mounts the current directory's `share` folder to `DATABASE_DIR` inside the container, creating a shared path between the host and the container.
##### 2. Testing and Quick Start
###### Testing
```
ccrr -mode test
```
This will test the necessary environment required by the process with an accompanying small sample data, which may take up to half an hour.
###### Quick Start
```
nohup ccrr -mode default \
        -normal [normal.bam] --normal-id [normal-id] \
        -tumor [tumor.bam] --tumor-id [tumor-id] \
        --genome-version hg19 -reference [hg19.fa] \
        -gender male -threads 30 -g 200 >log 2>&1 &
```
This will run the entire process in default mode, allowing for up to 30 threads for possible multi-threaded tasks, and a memory cap of 200GB. Processing a pair of tumor/normal control BAM files, each sized around 106GB, approximately takes 80 hours.

#### 3. Usage
##### Display Help
```
ccrr --help
```
##### Required Parameters
###### Mode
* fast: Only uses delly to obtain SV and CN data.
* default: Uses all included tools to obtain SV and CN data, then performs merging. SV results supported by at least two SV callers are retained by default.
* custom: Allows for the selection of tools, parameters, and merge method to be used.
* test: Tests the process with the accompanying small sample data, which may take about half an hour.
* clear: Clear the working directory for the task ID specified by the `-prefix` option. Please clear it before starting a new task.
```
  -mode {fast,custom,default,test,clear}        choose mode to run
```
###### Input and Information
Your input should be a pair of normal/tumor control whole-genome sequencing BAM files and their reference genome.  Supported reference genome versions are hg19 and hg38.
```
  -prefix                       task id  

  -normal NORMAL                normal bam
  --normal-id NORMAL_ID

  -tumor TUMOR                  tumor bam
  --tumor-id TUMOR_ID

  --genome-version {hg19,hg38}  Set the reference, hg19 or hg38
  -reference REFERENCE          reference fq
  -gender GENDER                male or female
```
##### Optional Parameters
###### Configuring Multithreading and Available Memory
If not set, the default memory allocation is 8GB, which may not suffice for the memory demands of certain steps. We recommend setting it higher.
Please note that the default number of threads is 8. Some software may not support multithreading acceleration, and for others, there might be a soft cap on the number of threads that can be effectively utilized, meaning that setting a higher number of threads may not result in the expected speed-up.
```
  -threads THREADS      Set the number of processes if possible
  -g G                  set the amount of available RAM If possible
```
###### Selecting Required Tools
In custom mode, you can freely choose which software to use for generating SV and CN data. The built-in software includes:
```
  -sequenza             use sequenza for cellularity and ploidy

  -delly                use delly for sv and cn 
  -manta                use manta for sv
  -svaba                use svaba for sv
  -gridss               use gridss for sv
  -lumpy                use lumpy for sv
  -soreca               use soreca for sv

  -sclust               use sclust for cn
  -purple               use purple for cn
  -cnvkit               use cnvkit for cn
```
###### Setting Quality Filtering
You can conveniently filter the results of each software based on quality before merging:
```
  --manta-filter MANTA_FILTER               Filter for manta
  --delly-filter DELLY_FILTER               Filter for delly sv
  --delly-cnvsize DELLY_CNVSIZE             min cnv size for delly
  --svaba-filter SVABA_FILTER               Filter for svaba
  --gridss-filter GRIDSS_FILTER             Filter for gridss
  --lumpy-filter LUMPY_FILTER               Filter for lumpy
```
###### Merging Method
When results from two different SV callers are adjacent in the genome and the distance between them is less than a specified threshold, they will be considered the same SV. The default threshold is 150bp.
```
  --sv-threshold SV_THRESHOLD
```
Select the method for merging results from different SV callers. If you wish to retain only the results supported by all the SV callers used, choose `intersection`; if you prefer to keep all results from all SV callers without duplicates, choose `union`; if you want to customize to retain results supported by X or more software tools, select `x-or-more`, and specify the number in `--sv-x`. If `X` is not specified, the default will be 3, meaning that results supported by two or more software tools will be retained.

```
  --sv-merge-method {intersection,union,x-or-more}
                        Choose a sv merging method: 
                        1. 'intersection': Merges only the SVs that are identified by all SV callers. 
                        2.'union': Merges all SVs identified by any of the SV callers. 
                        2. 'x-or-more': Merges SVs that are identified by at least x SV callers. if only one svcaller is prepared , then this parameter is irrelevant  
  --sv-x {1,2,3,4,5,6}  
                        Specify the x. This argument is required when '--merge-method' is set to 'x-or-more'. Must be among the provided input files. default=3
```
If you wish to prioritize a specific SV caller, setting `--sv-primary-caller` will retain all results outputted by it.
```
  --sv-primary-caller {manta,delly,svaba,gridss,lumpy,soreca}
                        Specify the primary SV caller to keep all of its result.
```
Setting `--cn-threshold` allows you to adjust the threshold, defining the maximum allowable distance for determining overlap among copy number change regions from different tools when merging copy number variant analysis results. The default threshold is 5000bp.
```
  --cn-threshold CN_THRESHOLD
                        threshold for determining cn, defaults to 5000bp
```
###### Complex Rearrangement Analysis
```
-complex COMPLEX      complex rearrangement analysis
```
Complex rearrangement analysis is conducted by default. If you only wish to obtain merged results, you can use `-complex False`.
##### Output, Rerunning, and History
`${WORK_DIR}/[task id]`will serve as the working directory, retaining the output results of each part. A summary of the complex rearrangement analysis can be found in `${WORK_DIR}/[task id]/complex/summary`.
Once a module is completed, it will be recorded in the `${WORK_DIR}/[task id]/history` file. If the process is unexpectedly interrupted, rerunning the entire process will skip the parts that have been successfully executed according to the records in the history file, resuming from the point of interruption.
Of course, you can manually modify this file to skip any steps you wish to bypass.
#### 4. Custom Execution
You can run each module step by step according to your analytical needs. For example:
##### Using `svmerge.py` to Merge SV Data
You can specify the output results from each SV caller as input files for merging.
```
  -manta MANTA          manta vcf result
  -delly DELLY          delly vcf result
  -svaba SVABA          svaba vcf result
  -gridss GRIDSS        GRIDSS vcf result
  -lumpy LUMPY          LUMPY vcf result
  -soreca SORECA        soreca result
```
Filter the output results of each SV caller based on quality.
```
  --manta-filter MANTA_FILTER
                        Filter threshold for manta
  --delly-filter DELLY_FILTER
                        Filter threshold for delly
  --svaba-filter SVABA_FILTER
                        Filter threshold for svaba
  --gridss-filter GRIDSS_FILTER
                        Filter threshold for GRIDSS
  --lumpy-filter LUMPY_FILTER
                        Filter threshold for LUMPY
```
Determine thresholds, merging methods, and specify a trusted SV caller as described previously.
```
  --threshold THRESHOLD
                        threshold for determination, defaults to 100bp

  --merge-method {intersection,union,x-or-more}
                        Choose a merging method: 1. 'intersection': Merges only the SVs that are identified by all SV callers. 
                        2. 'union': Merges all SVs identified by any of the SV callers. 
                        3.'x-or-more': Merges SVs that are identified by at least x SV callers. if only one svcaller is prepared , then this parameter is irrelevant

  --primary-caller {None,manta,delly,svaba,gridss,lumpy,soreca}
                        Specify the primary SV caller to keep all of its result.

  -x {1,2,3,4,5,6}      Specify the x. This argument is required when '--merge-method' is set to 'x-or-more'.
                        Must be among the provided input files.
```
Set the output path and enable multi-process execution.
```
  -o O                  output path
  -t T                  Set the number of processes
```
##### Use `consensus_cn.py` to merge CN data.
```
python ${SCRIPT_DIR}/consensus_cn.py  \
    -sclust SCLUST -delly DELLY -purple PURPLE -cnvkit CNVKIT \
    -ref hg19 -gender male \
    -o OUT 
```
Parameters
```
  -sclust SCLUST        sclust cn result
  -delly DELLY          delly cn result
  -purple PURPLE        purple cn result
  -cnvkit CNVKIT        cnvkit cn result
  --threshold THRESHOLD
                        threshold for determination, defaults to 5000bp

  -o O                  output path
  -ref REF              hg19 or hg38
```
##### Use `complex.py` to analyze complex rearrangements.
```
python ${SCRIPT_DIR}/complex.py  -prefix task_id \
        --tumor-id example -sv SV -cn CN \
        --genome-version hg19 -gender male \
        -shatterseek -starfish -gGnome -SA  -ctlpscanner \
        -threads 30 -g 200
```
Required inputs, the format of SV and CN files as shown in the examples.
```
http://life-bioinfo.tpddns.cn:13362/static/examplefile/custom_sv.bed
http://life-bioinfo.tpddns.cn:13362/static/examplefile/custom_cn.bed
```

```
  -prefix               task id
  --tumor-id TUMOR_ID
  -sv SV                sv input
  -cn CN                cn input
  --genome-version GENOME_VERSION
                        Set the reference, hg19 or hg38
  -gender {male,female}
                        gender
```
Select the tools for complex rearrangement analysis.
```
  -shatterseek          use shatterseek
  -starfish             use starfish
  -gGnome               use jabba and gGnome
  -SA                   use Seismic Amplification
  -ctlpscanner          use CTLPscanner
```
The AmpliconArchitect requires an input of BAM files.
```
  -AA                           use Amplicon Architect

  -normal NORMAL                normal bam
  --normal-id NORMAL_ID

  -tumor TUMOR                  tumor bam
  --tumor-id TUMOR_ID

```
Set the available memory and number of threads.
```
  -threads THREADS      Set the number of processes if possible
  -g G                  set the amount of available RAM If possible
```
