#!/opt/conda/envs/main/bin/python
import argparse
import json
import os
import re
import shutil
import subprocess
import config

def cmd(command, shell=False):
    if not shell:
        command = command.split()
    ret = subprocess.run(command, shell=shell)
    return ret.returncode == 0


def validate_file(filepath):
    if not os.path.exists(filepath):
        raise argparse.ArgumentTypeError(f"The file '{filepath}' does not exist!")
    else:
        return os.path.abspath(filepath)
    
def validate_id(id_str):
    if re.match(r'^\w+$', id_str):
        return id_str
    else:
        raise argparse.ArgumentTypeError("illegal id")

def rename_index_file(original_bam, new_bam_base):
    index_file1 = original_bam + '.bai'
    index_file2 = original_bam[:-4] + '.bai'

    if os.path.exists(index_file1):
        new_index_path1 = new_bam_base[:-4] + '.bai'
        new_index_path2 = new_bam_base + '.bai'
        shutil.move(index_file1, new_index_path1)
        shutil.copy(new_index_path1, new_index_path2)
    elif os.path.exists(index_file2):
        new_index_path1 = new_bam_base[:-4] + '.bai'
        new_index_path2 = new_bam_base + '.bai'
        shutil.move(index_file2, new_index_path1)
        shutil.copy(new_index_path1, new_index_path2)
    else:
        print("Index file not found. Please ensure the index file is present.")

def getarg():
    parser = argparse.ArgumentParser(description="main pipeline")

    parser.add_argument("-mode", type=str, default="default", choices=["fast", "custom", "default", "test", "clear"],
                        help="choose mode to run")
    parser.add_argument("-prefix", type=str, default="test")
    parser.add_argument("-normal", type=validate_file, help="normal bam")
    parser.add_argument("--normal-id", type=validate_id)
    parser.add_argument("-tumor", type=validate_file, help="tumor bam")
    parser.add_argument("--tumor-id", type=validate_id)
    
    parser.add_argument('-sequenza', type=bool, default=True,help="use sequenza for cellularity, ploidy and cn")
    parser.add_argument("-manta", action="store_true", help="use manta for sv")
    parser.add_argument("-delly", action="store_true", help="use delly for sv and cn")
    parser.add_argument("-svaba", action="store_true", help="use svaba for sv")
    parser.add_argument("-gridss", action="store_true", help="use gridss for sv")
    parser.add_argument("-lumpy", action="store_true", help="use lumpy for sv")
    parser.add_argument("-soreca", action="store_true", help="use soreca for sv")
    parser.add_argument("-purple", action="store_true", help="use purple for cn")
    parser.add_argument("-sclust", action="store_true", help="use sclust for cn")
    parser.add_argument("-cnvkit", action="store_true", help="use cnvkit for cn")

    parser.add_argument("--manta-filter", type=int, default=None, help="Filter for manta")
    parser.add_argument("--delly-filter", type=int, default=None, help="Filter for delly sv")
    parser.add_argument("--delly-cnvsize", type=int, default=1000, help="min cnv size for delly")
    parser.add_argument("--svaba-filter", type=int, default=None, help="Filter for svaba")
    parser.add_argument("--gridss-filter", type=int, default=None, help="Filter for gridss")
    parser.add_argument("--lumpy-filter", type=int, default=None, help="Filter for lumpy")
    parser.add_argument("--varscan-p", type=float, default=0.01, help="varscan p-value")
    parser.add_argument("--sclust-alpha", type=float, default=0.15, help="sclust alpha")
    parser.add_argument("--cellularity-ploidy-tool", type=str,choices=["sequenza", "purple"],default='sequenza')
    
    parser.add_argument("--genome-version", type=str, choices=["hg19", "hg38"],
                        help="Set the reference, hg19 or hg38")
    parser.add_argument("-reference", type=validate_file, help="reference fq")
    parser.add_argument("-gender", type=str, help="male or female")


    parser.add_argument("--sv-threshold", type=int, default=150,
                        help="threshold for determining breakpoints, defaults to 150bp")
    parser.add_argument("--sv-merge-method", type=str, choices=["intersection", "union", "x-or-more"],
                        default="x-or-more",
                        help="Choose a sv merging method: \n"
                             "1. 'intersection': Merges only the SVs that are identified by all SV callers.\n"
                             "2. 'union': Merges all SVs identified by any of the SV callers.\n"
                             "3. 'x-or-more': Merges SVs that are identified by at least x SV callers.\n"
                             "if only one svcaller is prepared , then this parameter is irrelevant")
    parser.add_argument("--sv-primary-caller", type=str,
                        choices=["manta", "delly", "svaba", "gridss", "lumpy", "soreca"],
                        help="Specify the primary SV caller to keep all of its result.")
    parser.add_argument("--sv-x", type=int, choices=[1, 2, 3, 4, 5, 6], default=3,
                        help="Specify the x. This argument is required when '--merge-method' is set to "
                             "'x-or-more'. Must be among the provided input files. default=3")
    parser.add_argument("--cn-threshold", type=int, default=5000,
                        help="threshold for determining cn, defaults to 5000bp")
    parser.add_argument("-threads", type=int, default=8, help="Set the number of processes if possible")
    parser.add_argument("-g", type=int, default=8, help="set the amount of available RAM If possible")

    parser.add_argument("-complex", type=bool, default=True, help="complex rearrangement analysis")
    args = parser.parse_args()

    if args.mode !="custom":
        sv_alltools = [tool for tool in ["manta", "delly", "svaba", "gridss", "lumpy", "soreca"] if tool in config.TOOLS_LIST]
        cn_alltools = [tool for tool in ["delly", "sclust", "cnvkit", "purple"] if tool in config.TOOLS_LIST]
        if args.mode == "fast":
            sv_tools = ["delly"]
            cn_tools = ["delly"]
            return sv_tools, cn_tools, args
        elif args.mode == "clear":
            return [], [], args
        elif args.mode == "default":
            sv_tools = sv_alltools
            cn_tools = cn_alltools
            return sv_tools, cn_tools, args
        elif args.mode == "test":
            sv_tools = sv_alltools
            cn_tools = cn_alltools
            args.normal = os.path.join(config.DATABASE_DIR,"data/test/normal/SRR5651359.bam")
            args.normal_id = "SRR5651359"
            args.tumor = os.path.join(config.DATABASE_DIR,"data/test/tumor/SRR5651358.bam")
            args.tumor_id = "SRR5651358"
            args.genome_version = "hg19"
            args.reference = os.path.join(config.DATABASE_DIR,"data/ref/hg19/hg19.genomic.fa")
            args.gender = "male"
        return sv_tools, cn_tools, args
    else:
        sv_tools = []
        cn_tools = []
        if args.manta:
            sv_tools.append("manta")
        if args.delly:
            sv_tools.append("delly")
            cn_tools.append("delly")
        if args.svaba:
            sv_tools.append("svaba")
        if args.gridss:
            sv_tools.append("gridss")
        if args.lumpy:
            sv_tools.append("lumpy")
        if args.soreca:
            sv_tools.append("soreca")
        if args.sclust:
            cn_tools.append("sclust")
        if args.cnvkit:
            cn_tools.append("cnvkit")
        if args.purple:
            cn_tools.append("purple")
            
        if len(sv_tools) == 0 or len(cn_tools) == 0:
            parser.error("you must select one sv caller and one cn caller at least")
            raise SystemExit

        if args.sv_merge_method == "x-or-more" and args.sv_x is None:
            parser.error("--sv-x is required when '--sv-merge-method' is set to 'x-or-more'")
            raise SystemExit

        if args.sv_merge_method == "x-or-more" and len(sv_tools) < args.sv_x:
            parser.error("x must be greater than the number of tools you have selected")
            raise SystemExit

        return sv_tools, cn_tools, args


def update_history(history_file, entry):
    with open(history_file, "a") as f:
        f.write(entry + "\n")


def main():
    sv_tools, cn_tools, args = getarg()
    workdir= os.path.join(config.WORK_DIR,args.prefix)
    history_file = os.path.join(workdir,"history")
    if args.mode == "clear":
        try:
            cmd("rm -rf "+workdir)
        except:
            pass
        return
    if not os.path.exists(workdir):
        cmd("mkdir -p "+ workdir)
    history = []
    cellularity = ploidy ='null'
    if os.path.exists(history_file):
        with open(history_file, "r") as f:
            history = f.read().splitlines()
    if args.normal and args.normal_id:
        original_dir = os.path.dirname(args.normal)
        new_normal_name = f"{args.normal_id}.bam"
        new_normal_path = os.path.join(original_dir, new_normal_name)
        shutil.move(args.normal, new_normal_path) 
        print(f"Normal file has been renamed to: {new_normal_path}")
        rename_index_file(args.normal, new_normal_path)
        args.normal = new_normal_path
    if args.tumor and args.tumor_id:
        original_dir = os.path.dirname(args.tumor)
        new_tumor_name = f"{args.tumor_id}.bam"
        new_tumor_path = os.path.join(original_dir, new_tumor_name)
        shutil.move(args.tumor, new_tumor_path)
        print(f"Tumor file has been renamed to: {new_tumor_path}")
        rename_index_file(args.tumor, new_tumor_path)
        args.tumor = new_tumor_path
    if 'sequenza' in config.TOOLS_LIST and 'sequenza' not in history and args.sequenza :
        if cmd("conda run --no-capture-output -n main bash "+os.path.join(config.SCRIPT_DIR,"sequenza.sh") + " -n "  + args.normal + " -t " + args.tumor + " -r " +args.reference + " -i " + str(args.prefix) + " -d " + str(args.threads)):
            history.append('sequenza')
            update_history(history_file, 'sequenza')
        else:
            cmd("rm -rf " + os.path.join(workdir, "sequenza"))
            return
    if "sequenza" in config.TOOLS_LIST and 'sequenza' in history and args.sequenza and args.cellularity_ploidy_tool == 'sequenza' :
        with open(os.path.join(workdir, "sequenza","cellularity_ploidy.json"), 'r') as file:
            data = json.load(file)
            cellularity = str(data['cellularity'][0])
            ploidy = str(data['ploidy'][0])
    if "delly" in sv_tools and "delly" not in history:
        if cmd("conda run --no-capture-output -n main bash "+os.path.join(config.SCRIPT_DIR,"delly.sh") + " -n "  + args.normal + " -t " + args.tumor + " -r " +
               args.reference + " -b " + args.genome_version + " -1 " + args.tumor_id + " -2 " + args.normal_id + " -z " + 
               str(args.delly_cnvsize) + " -x 2 -s 0.100000001 -i " + str(args.prefix)):
            history.append("delly")
            update_history(history_file, "delly")
        else:
            cmd("rm -rf "+ os.path.join(workdir,"delly"))
            print("delly error")
            return
    if "manta" in sv_tools and "manta" not in history:
        if cmd("conda run --no-capture-output -n py2 bash "+os.path.join(config.SCRIPT_DIR,"manta.sh") + " -n " + args.normal + " -t " + args.tumor + " -r " +
               args.reference + " -p " + str(args.threads) + " -g " + str(args.g) + " -i " + str(args.prefix)):
            history.append("manta")
            update_history(history_file, "manta")
        else:
            cmd("rm -rf "+ os.path.join(workdir,"manta"))
            print("manta error")
            return
    
    if "soreca" in sv_tools and "soreca" not in history:
        if cmd("conda run --no-capture-output -n main bash "+os.path.join(config.SCRIPT_DIR,"soreca.sh") + " -n " + args.normal + " -t " + args.tumor + " -r " + args.genome_version
               +" -i " + str(args.prefix)):
            history.append("soreca")
            update_history(history_file, "soreca")
        else:
            cmd("rm -rf "+ os.path.join(workdir,"soreca"))
            print("soreca error")
            return
        
    if "lumpy" in sv_tools and "lumpy" not in history:
        if cmd("conda run --no-capture-output -n py2 bash "+os.path.join(config.SCRIPT_DIR,"lumpy.sh")+" -n " + args.normal + " -t " + args.tumor+" -i " + str(args.prefix)):
            history.append("lumpy")
            update_history(history_file, "lumpy")
        else:
            cmd("rm -rf "+ os.path.join(workdir,"lumpy"))
            print("lumpy error")
            return
        
    if "svaba" in sv_tools and "svaba" not in history:
        if cmd("conda run --no-capture-output -n main bash "+os.path.join(config.SCRIPT_DIR,"svaba.sh")+" -n " + args.normal + " -t " + args.tumor + " -r " + args.reference
               + " -p " + str(args.threads)+" -i " + str(args.prefix) + " -b " + args.genome_version):
            history.append("svaba")
            update_history(history_file, "svaba")
        else:
            cmd("rm -rf "+ os.path.join(workdir,"svaba"))
            print("svaba error")
            return
        
    if "sclust" in cn_tools and "varscan" not in history:
        if cmd("conda run --no-capture-output -n main bash "+os.path.join(config.SCRIPT_DIR,"varscan.sh") + " -n " + args.normal + " -j " + str(
                args.g) + " -t " + args.tumor + " -i " + str(args.prefix) + " -b " +
               args.genome_version + " -r " + args.reference + " -p " + str(args.varscan_p)):
            history.append("varscan")
            update_history(history_file, "varscan")
        else:
            cmd("rm -rf " + os.path.join(workdir,"varscan"))
            print("varscan error")
            return
        
    if 'gridss' in sv_tools and 'gridss' not in history:
        if cmd("conda run --no-capture-output -n main bash "+os.path.join(config.SCRIPT_DIR,"gridss.sh")+ " -p " + str(args.threads) 
               + " -t " + args.tumor + " -r " + args.reference + " -i " + str(args.prefix) 
               + " -s True -b " + args.genome_version + " -n " + args.normal + " -1 " + args.normal_id + " -2 " + args.tumor_id + " -j " + str(args.g)):
            history.append('gridss')
            update_history(history_file, 'gridss')
        else:
            cmd("rm -rf " + os.path.join(workdir,"gridss"))
            print('gridss error')
            return
        
    if "sclust" in cn_tools and "varscan" in history and "sclust" not in history:
        if cmd("conda run --no-capture-output -n main bash "+os.path.join(config.SCRIPT_DIR,"sclust.sh")+" -a " + str(
                args.sclust_alpha)+" -i " + str(args.prefix) + " -n " + args.normal + " -e " + args.mode + " -t " + args.tumor
               + " -r " + args.genome_version):
            history.append("sclust")
            update_history(history_file, "sclust")
        else:
            cmd("rm -rf "+ os.path.join(workdir,"sclust"))
            print("sclust error")
            return
        
    if 'purple' in cn_tools and 'pre_purple' not in history:
        if cmd("conda run --no-capture-output -n main bash "+os.path.join(config.SCRIPT_DIR,"pre_purple.sh")+" -p " + str(
                args.threads) +" -i " + str(args.prefix) + " -r " + args.reference + " -b " + args.genome_version + " -n " +
               args.normal + " -1 " + args.normal_id + " -t " + args.tumor + " -2 " + args.tumor_id + " -j " + str(
            args.g)):
            history.append('pre_purple')
            update_history(history_file, 'pre_purple')
        else:
            cmd("rm -rf "+ os.path.join(workdir,"pre_purple"))
            print('pre_purple error')
            return
        
    if 'purple' in cn_tools and "pre_purple" in history and 'purple' not in history:
        if 'gridss' in sv_tools and 'gridss' in history:
            if cmd("conda run --no-capture-output -n main bash "+os.path.join(config.SCRIPT_DIR,"purple.sh")+ " -p " + str(
                    args.threads) + " -r " + args.reference + " -b " + args.genome_version + " -n " +
                args.normal + " -1 " + args.normal_id + " -t " + args.tumor + " -2 " + args.tumor_id + " -j " + str(
                args.g)+" -i " + str(args.prefix) + " -g True"):
                history.append("purple")
                update_history(history_file, 'purple')
            else:
                cmd("rm -rf "+ os.path.join(workdir,"purple"))
                print('purple error')
                return
        else:
            if cmd("conda run --no-capture-output -n main bash "+os.path.join(config.SCRIPT_DIR,"purple.sh")+ " -p " + str(
                    args.threads) + " -r " + args.reference + " -b " + args.genome_version + " -n " +
                args.normal + " -1 " + args.normal_id + " -t " + args.tumor + " -2 " + args.tumor_id + " -j " + str(
                args.g)+" -i " + str(args.prefix) + " -g False"):
                history.append("purple")
                update_history(history_file, 'purple')
            else:
                cmd("rm -rf "+ os.path.join(workdir,"purple"))
                print('purple error')
                return
    if "purple" in config.TOOLS_LIST and 'purple' in history and args.purple and args.cellularity_ploidy_tool =='purple' :
        with open(os.path.join(workdir, "sequenza","output",args.tumor_id+".purple.purity.tsv"), 'r') as file:
            headers = file.readline().strip().split('\t')
            values = file.readline().strip().split('\t')
            cellularity, ploidy = values[headers.index('purity')], values[headers.index('ploidy')]

    if "cnvkit" in cn_tools and "cnvkit" not in history:
        if cmd("conda run --no-capture-output -n main bash "+os.path.join(config.SCRIPT_DIR,"cnvkit.sh") +
               " -i " + str(args.prefix) + " -p " + str(args.threads)  + " -n " + args.normal + " -r " + args.reference + " -t " + args.tumor
               + " -b " + str(args.genome_version) + " -u " + str(cellularity) + " -l " + str(ploidy)) :
            history.append("cnvkit")
            update_history(history_file, "cnvkit")
        else:
            cmd("rm -rf "+ os.path.join(workdir,"cnvkit"))
            print("cnvkit error")
            return
    
    all_sv_tools_ready = all(tool in history for tool in sv_tools)
    all_cn_tools_ready = all(tool in history for tool in cn_tools)

    if not all_sv_tools_ready or not all_cn_tools_ready:
        print("Some results are not ready")
        return
    if "svmerge" not in history:
        cmd("rm -rf "+os.path.join(workdir,"svmerge"))
        cmd("mkdir "+os.path.join(workdir,"svmerge"))
        tools_cmd = ""
        if "delly" in sv_tools and cmd("cp "+ os.path.join(workdir,"delly","delly.sv.somatic.pre.vcf") + " "+os.path.join(workdir,"svmerge")):
            tools_cmd = tools_cmd + " -delly "+os.path.join(workdir,"svmerge","delly.sv.somatic.pre.vcf")
        if 'manta' in sv_tools:
            if cmd("gzip -dk "+ os.path.join(workdir,"manta","results","variants","somaticSV.vcf.gz")) and \
                    cmd("cp "+ os.path.join(workdir,"manta","results","variants","somaticSV.vcf") + " " + os.path.join(workdir,"svmerge")):
                tools_cmd = tools_cmd + " -manta "+os.path.join(workdir,"svmerge","somaticSV.vcf")
        if "svaba" in sv_tools: 
            if cmd("cp " + os.path.join(workdir, "svaba", "svaba.svaba.somatic.sv.vcf") + " " + os.path.join(workdir, "svmerge")):
                tools_cmd = tools_cmd + " -svaba " + os.path.join(workdir, "svmerge", "svaba.svaba.somatic.sv.vcf")
        if "gridss" in sv_tools:
            if cmd("gzip -dk " + os.path.join(workdir, "gridss", "gripss_output", args.tumor_id + ".gripss.filtered.vcf")) and \
                    cmd("cp " + os.path.join(workdir, "gridss", "gripss_output", args.tumor_id + ".gripss.filtered.vcf") + " " + os.path.join(workdir, "svmerge")):
                tools_cmd = tools_cmd + " -gridss " + os.path.join(workdir, "svmerge", args.tumor_id + ".gripss.filtered.vcf")
        if "lumpy" in sv_tools and cmd("cp " + os.path.join(workdir, "lumpy", "lumpy.gt.vcf") + " " + os.path.join(workdir, "svmerge")):
            tools_cmd = tools_cmd + " -lumpy " + os.path.join(workdir, "svmerge", "lumpy.gt.vcf")
        if "soreca" in sv_tools and cmd("cp " + os.path.join(workdir, "soreca", "soreca_unsnarl.txt") + " " + os.path.join(workdir, "svmerge")):
            tools_cmd = tools_cmd + " -soreca " + os.path.join(workdir, "svmerge", "soreca_unsnarl.txt")

        if cmd("conda run --no-capture-output -n main python "+os.path.join(config.SCRIPT_DIR,"svmerge.py")+" --merge-method " + args.sv_merge_method + " -x "+ str(args.sv_x) + tools_cmd +
               " --primary-caller " + str(args.sv_primary_caller) + " -t " + str(args.threads) + " -o "+os.path.join(workdir,"svmerge")+" --threshold " + str(args.sv_threshold)):
            history.append("svmerge")
            cmd("conda run --no-capture-output -n main Rscript "+os.path.join(config.SCRIPT_DIR,"draw","ccrr_svupset.R") +" -i "+os.path.join(workdir,"svmerge","sv_merged.bed") + " -o " + os.path.join(workdir,"svmerge","sv_merged.pdf") )
            update_history(history_file, "svmerge")
        else:
            print("sv merge error")
            return
    
    if "cnmerge" not in history:
        cmd("rm -rf "+os.path.join(workdir,"cnmerge"))
        cmd("mkdir "+os.path.join(workdir,"cnmerge"))
        tools_cmd = " "
        if "delly" in cn_tools and cmd("cp "+ os.path.join(workdir,"delly","segmentation.bed") + " " + os.path.join(workdir,"cnmerge")):
            tools_cmd = tools_cmd + " -delly "+ os.path.join(workdir,"cnmerge","segmentation.bed")
        if "purple" in cn_tools and cmd("cp " + os.path.join(workdir, "purple", "output", str(args.tumor_id) + ".purple.cnv.somatic.tsv") + " " + os.path.join(workdir, "cnmerge")):
            tools_cmd = tools_cmd + " -purple " + os.path.join(workdir, "cnmerge", str(args.tumor_id) + ".purple.cnv.somatic.tsv")
        if "sclust" in cn_tools and cmd("cp " + os.path.join(workdir, "sclust", "sclust_final_iCN.seg") + " " + os.path.join(workdir, "cnmerge")):
            tools_cmd = tools_cmd + " -sclust " + os.path.join(workdir, "cnmerge", "sclust_final_iCN.seg")
        if "cnvkit" in cn_tools and cmd("cp " + os.path.join(workdir, "cnvkit", str(args.prefix)+".cnvkit.bed") + " " + os.path.join(workdir, "cnmerge")):
            tools_cmd = tools_cmd + " -cnvkit " + os.path.join(workdir, "cnmerge", str(args.prefix)+".cnvkit.bed")
        if "sequenza" in config.TOOLS_LIST and args.sequenza and  cmd("cp " + os.path.join(workdir, "sequenza", str(args.prefix)+"_segments.txt") + " " + os.path.join(workdir, "cnmerge")):
            tools_cmd = tools_cmd + " -sequenza " + os.path.join(workdir, "cnmerge", str(args.prefix)+"_segments.txt")
        if cmd("conda run --no-capture-output -n main python " + os.path.join(config.SCRIPT_DIR,"consensus_cn.py") + tools_cmd + " -o " + os.path.join(workdir,"cnmerge") + " -ref " + args.genome_version):
            history.append("cnmerge")
            update_history(history_file, "cnmerge")
            cmd("conda run --no-capture-output -n main Rscript "+os.path.join(config.SCRIPT_DIR,"draw","ccrr_cn_bias_vio.R") +" -i "+os.path.join(workdir,"cnmerge","consensus_cn_statistic.bed") + " -o " + os.path.join(workdir,"cnmerge"))
        else:
            print("cn merge error")
            return

    if args.complex and ("complex" not in history or "draw" not in history):
        if "complex" not in history:
            cmd("rm -rf " + os.path.join(workdir,"complex"))
            cmd("mkdir " + os.path.join(workdir,"complex"))
        svpath = os.path.join(workdir,"svmerge","sv_merged.bed")
        cnpath = os.path.join(workdir,"cnmerge","consensus_cn.bed")
        
        if args.mode == "test":
            svpath = os.path.join(config.DATABASE_DIR,"data/test/custom_sv.bed")
            cnpath = os.path.join(config.DATABASE_DIR,"data/test/custom_cn.bed")

        if cmd("conda run --no-capture-out -n main python " + os.path.join(config.SCRIPT_DIR,"complex.py") +" -prefix " + str(args.prefix) + " --tumor-id " + args.tumor_id + " -sv " + svpath
               + " -cn " + cnpath + " --genome-version " + args.genome_version + " -tumor " +args.tumor +" -normal "+args.normal + " -reference " + args.reference + " -gender " + args.gender +
               " -shatterseek -starfish -gGnome -SA -AA -ctlpscanner -threads " + str(args.threads) + " -g " + str(args.g)+ " -purity " + str(cellularity)+ " -ploidy " + str(ploidy)):
            history.append("complex")
            update_history(history_file, "complex")
        else:
            print("complex analysis error")
            return


if __name__ == "__main__":
    main()