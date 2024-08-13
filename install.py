import argparse
import os
import tarfile
import requests
import zipfile
def validate_file(filepath):
    if not os.path.exists(filepath):
        raise argparse.ArgumentTypeError(f"The file '{filepath}' does not exist!")
    else:
        return os.path.abspath(filepath)


def getarg():
    parser = argparse.ArgumentParser(description="Select the software to be installed, then generate a Dockerfile.")
    parser.add_argument('-sequenza', action='store_true', help="use sequenza for cellularity and ploidy")
    parser.add_argument("-manta", action='store_true', help="use manta for sv")
    parser.add_argument("-delly", action='store_true', help="use delly for sv and cn")
    parser.add_argument("-svaba", action='store_true', help="use svaba for sv")
    parser.add_argument("-gridss", action='store_true', help="use gridss for sv")
    parser.add_argument("-lumpy", action='store_true', help="use lumpy for sv")
    parser.add_argument("-soreca", action='store_true', help="use soreca for sv")
    parser.add_argument("-purple", action='store_true', help="use purple for cn")
    parser.add_argument("-sclust", action='store_true', help="use sclust for cn")
    parser.add_argument("-cnvkit", action='store_true', help="use cnvkit for cn")
    parser.add_argument("-ref", type=str,choices=['hg19',"hg19&hg38"],default='hg19')
    
    args = parser.parse_args()
    tools = []
    sv_tools = ['manta', 'delly', 'svaba', 'gridss', 'lumpy', 'soreca']
    cn_tools = ['delly', 'purple', 'sclust', 'cnvkit','sequenza']
    sv_selected = any(getattr(args, tool) for tool in sv_tools)
    cn_selected = any(getattr(args, tool) for tool in cn_tools)

    if not sv_selected or not cn_selected:
        parser.error("At least one SV tool and one CN tool must be selected.")
    for arg in vars(args):
        if getattr(args, arg) == True:
            tools.append(arg)

    return tools,args

def dockerfile(tools, args):
    tools_string = ",".join(tools)
    dockerfile_content = f"""FROM condaforge/mambaforge:latest
ARG GITHUB_PAT
ARG SCRIPT_DIR="/home/0.script"
ARG TOOL_DIR="/home/1.tools"
ARG DATABASE_DIR="/home/2.share"
ARG WORK_DIR="/home/3.wd"
ARG TOOLS_LIST="{tools_string}"

ENV SCRIPT_DIR=$SCRIPT_DIR TOOL_DIR=$TOOL_DIR DATABASE_DIR=$DATABASE_DIR WORK_DIR=$WORK_DIR TOOLS_LIST=$TOOLS_LIST

RUN mkdir -p ${{SCRIPT_DIR}} && \\
    echo "SCRIPT_DIR = \\"${{SCRIPT_DIR}}\\"" > ${{SCRIPT_DIR}}/config.py && \\
    echo "TOOL_DIR = \\"${{TOOL_DIR}}\\"" >> ${{SCRIPT_DIR}}/config.py && \\
    echo "DATABASE_DIR = \\"${{DATABASE_DIR}}\\"" >> ${{SCRIPT_DIR}}/config.py && \\
    echo "WORK_DIR = \\"${{WORK_DIR}}\\"" >> ${{SCRIPT_DIR}}/config.py && \\
    echo "TOOLS_LIST = \\"${{TOOLS_LIST}}\\"" >> ${{SCRIPT_DIR}}/config.py && \\
    conda create -n main && \\
    echo "source activate main" > ~/.bashrc && \\
    . /opt/conda/etc/profile.d/conda.sh && \\
    conda activate main && \\
    mamba install python=3.10.13 r-base=4.2.3 conda-forge::pyvcf=0.6.8 conda-forge::pandas=2.2.1 conda-forge::r-devtools=2.4.5 vim\\
"""
    if args.delly:
        dockerfile_content = dockerfile_content + """        bioconda::delly=1.1.8 compbiocore::bcftools=1.6 \\
"""
    if args.cnvkit:
        dockerfile_content = dockerfile_content + """        bioconda::cnvkit=0.9.10 \\
"""
    if args.gridss:
        dockerfile_content = dockerfile_content + """        bioconda::gridss=2.13.2 bioconda::hmftools-gripss=2.3.2 \\
"""
    if args.purple:
        dockerfile_content = dockerfile_content + """        bioconda::hmftools-cobalt=1.13=hdfd78af_1 bioconda::hmftools-amber=3.9=hdfd78af_1 bioconda::hmftools-sage=3.2.3=hdfd78af_0 bioconda::hmftools-purple=3.7.1=hdfd78af_0 \\
"""
    if args.sclust:
        dockerfile_content = dockerfile_content + """        bioconda::samtools=1.18 conda-forge::openjdk=20.0.2  bioconda::varscan=2.4.6 \\
"""
    if args.sequenza:
        dockerfile_content = dockerfile_content + """        bioconda::r-sequenza=3.0.0  r-jsonlite=1.8.8 bioconda::sequenza-utils=3.0.0 \
bioconda::bioconductor-genomeinfodbdata=1.2.9 bioconda::bioconductor-genomeinfodb=1.34.9 conda-forge::r-iotools=0.3_2 \\
"""
    if args.svaba:
        dockerfile_content = dockerfile_content + """        bioconda::svaba=1.1.0 bioconda::bwa=0.7.18 bioconda::samtools=1.18 \\
"""
    dockerfile_content = dockerfile_content + """        bioconda::bioconductor-biocgenerics=0.44.0 bioconda::bioconductor-graph=1.76.0 \\
        bioconda::bioconductor-s4vectors=0.36.0 bioconda::bioconductor-genomicranges=1.50.0 \\
        bioconda::bioconductor-iranges=2.32.0 conda-forge::r-gridextra=2.3 conda-forge::r-ggplot2=3.5.0 \\
        bioconda::bioconductor-consensusclusterplus=1.62.0 r-foreach=1.5.2 r-neuralnet=1.44.2 r-plyr=1.8.9 r-data.table=1.15.2 \\
        r-mass=7.3_60.0.1 r-gridextra=2.3 r-dplyr=1.1.4 r-factoextra=1.0.7 r-dendextend=1.17.1 r-gplots=3.1.3.1 r-ggpubr=0.6.0 \\
        r-reshape2=1.4.4 r-cowplot=1.1.3 r-patchwork=1.2.0 r-cairo=1.6_1 r-ggforce=0.4.2 r-testthat=3.2.1 \\
        conda-forge::r-markovchain=0.9.5 bioconda::bioconductor-copynumber=1.38.0  bioconda::bioconductor-rtracklayer=1.58.0 \\
        bioconda::ampliconsuite mosek::mosek=10.1.28 bioconda::bioconductor-complexheatmap=2.14.0 conda-forge::r-complexupset conda-forge::r-circlize=0.4.16 -y && \\
    conda create -n py2 && . /opt/conda/etc/profile.d/conda.sh && \\
        conda activate py2 && mamba install python=2.7.15 \\
"""
    if args.manta:
        dockerfile_content = dockerfile_content + """        bioconda::manta=1.6.0  \\
"""    
    if args.lumpy:
        dockerfile_content = dockerfile_content + """        compbiocore::pysam=0.12.0.1 compbiocore::samtools=1.3.1 bioconda::svtyper=0.6.1 bioconda::lumpy-sv=0.3.1 \\
"""
    dockerfile_content = dockerfile_content + "        -y\n"
    
    if args.delly:
        dockerfile_content = dockerfile_content + """COPY tools/delly ${TOOL_DIR}/delly
"""
    if args.gridss:
        dockerfile_content = dockerfile_content + """COPY tools/gridss ${TOOL_DIR}/gridss
"""
    if args.sclust:
        dockerfile_content = dockerfile_content + """COPY tools/sclust ${TOOL_DIR}/sclust
"""
    if args.soreca:
        dockerfile_content = dockerfile_content + """COPY tools/soreca ${TOOL_DIR}/soreca
"""
    dockerfile_content = dockerfile_content + """COPY tools/AmpliconSuite ${TOOL_DIR}/AmpliconSuite
COPY tools/ctlpscanner ${TOOL_DIR}/ctlpscanner
COPY tools/jabba ${TOOL_DIR}/jabba
COPY tools/SA ${TOOL_DIR}/SA
COPY tools/Starfish ${TOOL_DIR}/Starfish
COPY gurobi.lic /opt/gurobi/gurobi.lic
COPY mosek.lic /root/mosek/mosek.lic
COPY script ${SCRIPT_DIR}
RUN . /opt/conda/etc/profile.d/conda.sh && \\
    conda activate main && \\
    set -e && \\
    export GITHUB_PAT=${GITHUB_PAT} && \\
    R -e "options(timeout=360);devtools::install_github('parklab/ShatterSeek', upgrade = 'never') " && \\
    R -e "options(timeout=360);install.packages('${TOOL_DIR}/jabba/gurobi1003/linux64/R/gurobi_10.0-3_R_4.2.0.tar.gz', repos = NULL, type = 'source')" && \\
    R -e "options(timeout=360);install.packages('slam', repos='http://cran.rstudio.com/')" && \\
    R CMD INSTALL ${TOOL_DIR}/Starfish/ShatterSeeky_0.4.tar.gz && \\
    R -e "options(timeout=360);devtools::install_github('mskilab-org/JaBbA', upgrade = 'never')" && \\
    R -e "options(timeout=360);devtools::install_github('yanglab-computationalgenomics/Starfish',upgrade = 'never')" && \\
    R -e "options(timeout=360);devtools::install_local('${TOOL_DIR}/jabba/gUtils', dependencies = FALSE,force=TRUE)" && \\
    R -e "options(timeout=360);devtools::install_local('${TOOL_DIR}/jabba/gGnome', dependencies = FALSE,force=TRUE)"

RUN /bin/bash -c "source /opt/conda/etc/profile.d/conda.sh && \\
    conda activate main && \\
    source ${TOOL_DIR}/AmpliconSuite/install.sh --finalize_only --data_repo_loc ${DATABASE_DIR}/database/AA"

RUN ln -sf ${SCRIPT_DIR}/ccrr.py /usr/local/bin/ccrr && \\
    chmod +x ${SCRIPT_DIR}/ccrr.py

CMD /bin/bash -c " """
    if args.lumpy:
        dockerfile_content = dockerfile_content + """ln -sf /opt/conda/envs/py2/lib/libcrypto.so.1.1 /opt/conda/envs/py2/lib/libcrypto.so.1.0.0 && \\
"""
    if args.sclust:
        dockerfile_content = dockerfile_content + """      chmod +x ${TOOL_DIR}/sclust/bin/Sclust && ln -sf ${DATABASE_DIR}/database/annotation ${TOOL_DIR}/sclust/annotation && \\
"""
    if args.soreca:
        dockerfile_content = dockerfile_content + """       chmod +x ${TOOL_DIR}/soreca/bin/soreca && chmod +x ${TOOL_DIR}/soreca/src/blat/blat && ln -sf ${DATABASE_DIR}/database/annotation ${TOOL_DIR}/soreca/annotation && \\
"""
    dockerfile_content = dockerfile_content + """        . /opt/conda/etc/profile.d/conda.sh && conda activate main && tail -f /dev/null"
"""
    with open('Dockerfile', 'w') as file:
        file.write(dockerfile_content)
    with open('.dockerignore', 'w') as file:
        file.write("share")

def download_file(url, filename):
    with requests.get(url, stream=True) as response:
        response.raise_for_status()
        with open(filename, 'wb') as f:
            for chunk in response.iter_content(chunk_size=32768):
                if chunk:
                    f.write(chunk)
    print(filename + " downloaded successfully")

def unzip_tar_gz(tar_path, extract_to):
    if not os.path.exists(extract_to):
        os.makedirs(extract_to)
    with tarfile.open(tar_path, 'r:gz') as tar:
        tar.extractall(path=extract_to)
        print("Files successfully extracted to", extract_to)
    os.remove(tar_path)
    print("Original tar.gz file removed:", tar_path)

def unzip(zip_path, extract_to):
    if not os.path.exists(extract_to):
        os.makedirs(extract_to)
    with zipfile.ZipFile(zip_path, 'r') as zip_ref:
        zip_ref.extractall(extract_to)
        print("Files successfully extracted to", extract_to)
    os.remove(zip_path)
    print("Original zip file removed:", zip_path)

def database(args):
    if args.gridss or args.purple:
        gridss_database=os.path.join("share","database","gridss_database")
        os.makedirs(gridss_database, exist_ok=True)
        if 'hg19' in args.ref:
            filename = "gridss_37.zip"
            download_file("http://life-bioinfo.tpddns.cn:13362/download_file/gridss_37.zip", \
                os.path.join(gridss_database,filename))
            unzip(os.path.join(gridss_database,filename), os.path.join(gridss_database,"37"))
        if 'hg38' in args.ref:
            filename = "gridss_38.zip"
            download_file("http://life-bioinfo.tpddns.cn:13362/download_file/gridss_38.zip", \
                os.path.join(gridss_database,filename))
            unzip(os.path.join(gridss_database,filename), os.path.join(gridss_database,"38"))
    if args.sclust or args.soreca:
        annotation_database=os.path.join("share","database","annotation")
        os.makedirs(annotation_database, exist_ok=True)
        if 'hg19' in args.ref:
            filename = "annotation.hg19.zip"
            download_file("http://life-bioinfo.tpddns.cn:13362/download_file/annotation.hg19.zip", \
                os.path.join(annotation_database,filename))
            unzip(os.path.join(annotation_database,filename), annotation_database)
        if 'hg38' in args.ref:
            filename = "annotation.hg38.zip"
            download_file("http://life-bioinfo.tpddns.cn:13362/download_file/annotation.hg38.zip", \
                os.path.join(annotation_database,filename))
            unzip(os.path.join(annotation_database,filename), annotation_database)
    if args.delly:
        delly_database=os.path.join("share","database","delly_map")
        os.makedirs(delly_database, exist_ok=True)
        filename = "delly_map.zip"
        download_file("http://life-bioinfo.tpddns.cn:13362/download_file/delly_map.zip", \
            os.path.join(delly_database,filename))
        unzip(os.path.join(delly_database,filename),delly_database)
    if args.svaba:
        dbsnp_database=os.path.join("share","database","dbsnp")
        os.makedirs(dbsnp_database, exist_ok=True)
        filename = "dbsnp.zip"
        download_file("http://life-bioinfo.tpddns.cn:13362/download_file/dbsnp.zip", \
            os.path.join(dbsnp_database,filename))
        unzip(os.path.join(dbsnp_database,filename),dbsnp_database)

    aa_database=os.path.join("share","database","AA","data_repo")
    os.makedirs(aa_database, exist_ok=True)
    if 'hg19' in args.ref:
        filename = "hg19.tar.gz"
        download_file("https://datasets.genepattern.org/data/module_support_files/AmpliconArchitect/hg19.tar.gz", \
            os.path.join(aa_database,filename))
        unzip_tar_gz(os.path.join(aa_database,filename), aa_database)
    if 'hg38' in args.ref:
        filename = "GRCh38.tar.gz"
        download_file("https://datasets.genepattern.org/data/module_support_files/AmpliconArchitect/GRCh38.tar.gz", \
            os.path.join(aa_database,filename))
        unzip_tar_gz(os.path.join(aa_database,filename), aa_database)
        
    return True


def main():
    tools, args = getarg()
    if database(args):
        dockerfile(tools,args)


if __name__ == "__main__":
    main()