import argparse
import os
import shutil
import subprocess
import config


def cmd(command, shell=False):
    if not shell:
        command = command.split()
    ret = subprocess.run(command, shell=shell)
    return ret.returncode == 0


def update_history(history_file, entry):
    with open(history_file, 'a') as f:
        f.write(entry + "\n")


def validate_file(filepath):
    if not os.path.exists(filepath):
        raise argparse.ArgumentTypeError(f"The file '{filepath}' does not exist!")
    else:
        return os.path.abspath(filepath)


def getarg():
    parser = argparse.ArgumentParser(description="complex SV pipeline")
    parser.add_argument('-prefix', type=str, default='test')
    parser.add_argument('-normal', type=validate_file, help="normal bam")
    parser.add_argument('--normal-id', type=str)
    parser.add_argument('-tumor', type=validate_file, help="tumor bam")
    parser.add_argument('--tumor-id', type=str)
    parser.add_argument('-sv', type=validate_file, help='sv input')
    parser.add_argument('-cn', type=validate_file, help='cn input')
    parser.add_argument('--genome-version', type=str, help="Set the reference, hg19 or hg38")
    parser.add_argument('-reference', type=validate_file, help="reference fq")
    parser.add_argument('-gender', default='male', choices=["male", "female"], help="gender")
    parser.add_argument('-purity', type=str,default="null")
    parser.add_argument('-ploidy', type=str,default="null")
    parser.add_argument('-shatterseek', action='store_true', help="use shatterseek")
    parser.add_argument('-starfish', action='store_true', help="use starfish")
    parser.add_argument('-gGnome', action='store_true', help="use jabba and gGnome")
    parser.add_argument('-AA', action='store_true', help="use Amplicon Architect")
    parser.add_argument('-SA', action='store_true', help="use Seismic Amplification")
    parser.add_argument('-ctlpscanner', action='store_true', help="use CTLPscanner")
    
    parser.add_argument('-threads', type=int, default=8, help="Set the number of processes if possible")
    parser.add_argument('-g', type=int, default=8, help="set the amount of available RAM If possible")

    args = parser.parse_args()

    complex_tools = [tool for arg, tool in [
        (args.shatterseek, 'shatterseek'),
        (args.starfish, 'starfish'),
        (args.gGnome, 'gGnome'),
        (args.AA, 'AA'),
        (args.SA, 'SA'),
        (args.ctlpscanner, 'ctlpscanner')
    ] if arg]

    return complex_tools, args


def main():
    complex_tools, args = getarg()
    work_path= os.path.join(config.WORK_DIR,args.prefix)
    history = []
    history_file = os.path.join(work_path, "history")
    if not os.path.exists(work_path):
        cmd('mkdir ' + work_path)
    complex_path = os.path.join(work_path, 'complex')
    if not os.path.exists(complex_path):
        cmd('mkdir ' + complex_path)
    if os.path.exists(history_file):
        with open(history_file, 'r') as f:
            history = f.read().splitlines()
    if "shatterseek" in complex_tools and "shatterseek" not in history:
        if cmd("conda run --no-capture-output -n main Rscript "+os.path.join(config.SCRIPT_DIR,"shatterseek.R") +" "+args.sv + " " + args.cn + " " + args.genome_version + " " + os.path.join(
                complex_path, "shatterseek")):
            history.append('shatterseek')
            update_history(history_file, 'shatterseek')
    if 'starfish' in complex_tools and 'starfish' not in history:
        if cmd("conda run --no-capture-output -n main Rscript "+os.path.join(config.SCRIPT_DIR,"starfish.R")+" "+ args.sv + " " + args.cn + " " + args.genome_version + " " + args.tumor_id + " " + args.gender + " " + os.path.join(
                complex_path, "starfish")):
            history.append('starfish')
            update_history(history_file, 'starfish')

    if 'gGnome' in complex_tools and 'gGnome' not in history:
        if cmd("conda run --no-capture-output -n main bash "+os.path.join(config.SCRIPT_DIR,"gGnome.sh")+ " -p " + str(
                    args.threads) + " -j " + str(args.g) + ' -s ' + args.sv + ' -c ' + args.cn + ' -o ' + os.path.join(
                complex_path, "jabba")+ " -u " + str(args.purity) + " -l " + str(args.ploidy)):
                history.append('gGnome')
                update_history(history_file, 'gGnome')
        cmd("rm -rf " + os.path.join(complex_path, "jabba"))

        
    if 'ctlpscanner' in complex_tools and 'ctlpscanner' not in history:
        if cmd("conda run --no-capture-output -n main bash "+os.path.join(config.SCRIPT_DIR,"ctlpscanner.sh")+ " -c " + args.cn + " -b " + args.genome_version + " -1 " + args.tumor_id + " -o " + os.path.join(
                complex_path, "ctlpscanner")):
            history.append('ctlpscanner')
            update_history(history_file, 'ctlpscanner')

    if 'AA' in complex_tools and 'AA' not in history:
        if cmd("conda run --no-capture-output -n main bash "+os.path.join(config.SCRIPT_DIR,"AA.sh") + " -n " + args.normal + " -t " + args.tumor + " -r " + args.genome_version + " -d " + str(
                    args.threads) + " -p " + args.tumor_id + " -c " +args.cn + " -i "+args.prefix):
            history.append('AA')
            update_history(history_file, 'AA')
        else:
            history.append('AA')
            update_history(history_file, 'AA')
            return

    if 'SA' in complex_tools and 'SA' not in history:
        if cmd("conda run --no-capture-output -n main Rscript "+ os.path.join(config.SCRIPT_DIR,"seismic_amplification.R") + " "+ args.cn + " " + args.sv + " " + args.genome_version + " " + os.path.join(
                complex_path, "SA")+" "+config.TOOL_DIR):
            history.append('SA')
            update_history(history_file, 'SA')


    if 'draw' not in history:
        complex_blank_path =  os.path.join(config.DATABASE_DIR,"database/complex_blank/")

        def copy_file(src, dst):
            os.makedirs(os.path.dirname(dst), exist_ok=True)
            shutil.copy(src, dst)

        shatterseek_file = os.path.join(complex_path, "shatterseek", "chromothripsis_summary.csv")
        if not os.path.exists(shatterseek_file):
            copy_file(os.path.join(complex_blank_path, "shatterseek", "chromothripsis_summary.csv"), shatterseek_file)

        ctlpscanner_file = os.path.join(complex_path, "ctlpscanner", "CTLPRegion.txt")
        if not os.path.exists(ctlpscanner_file):
            copy_file(os.path.join(complex_blank_path, "ctlpscanner", "CTLPRegion.txt"), ctlpscanner_file)

        sa_file = os.path.join(complex_path, "SA", "SA_amplicons.csv")
        if not os.path.exists(sa_file):
            copy_file(os.path.join(complex_blank_path, "SA", "SA_amplicons.csv"), sa_file)

        gGnome_file = os.path.join(complex_path, "gGnome", "event_footprints.txt")
        if not os.path.exists(gGnome_file):
            copy_file(os.path.join(complex_blank_path, "gGnome", "event_footprints.txt"), gGnome_file)
        aa_file = os.path.join(complex_path, "AA", args.tumor_id + "_classification",args.tumor_id + "_result_data.json")
        if not os.path.exists(gGnome_file):
            copy_file(os.path.join(complex_blank_path, "AA", "result_data.json"), aa_file)

        starfish_file = os.path.join(complex_path, "starfish", args.tumor_id + "_connected_CGR_event.csv")
        if not os.path.exists(starfish_file):
            copy_file(os.path.join(complex_blank_path, "starfish", "tumor_connected_CGR_event.csv"), starfish_file)

        if args.genome_version == 'hg19' and args.gender == 'male' and cmd(
                "conda run --no-capture-output -n main Rscript "+ os.path.join(config.SCRIPT_DIR,"draw","ccrr_draw.R")+" --ref "+ os.path.join(config.SCRIPT_DIR,"draw/hg19_man.txt")+" --cn "
                + args.cn + " --shatterseek " + os.path.join(complex_path, "shatterseek", "chromothripsis_summary.csv")
                + " --region " + os.path.join(config.SCRIPT_DIR,"draw/region/hg19_man.region")
                + " --ctlpscanner " + os.path.join(complex_path, "ctlpscanner", "CTLPRegion.txt")
                + " --sa " + os.path.join(complex_path, "SA", "SA_amplicons.csv")
                + " --aa " + os.path.join(complex_path, "AA", args.tumor_id + "_classification",args.tumor_id + "_result_data.json")
                + " --gGnome " + os.path.join(complex_path, "gGnome", "event_footprints.txt")
                + " --starfish " + os.path.join(complex_path, "starfish", args.tumor_id + "_connected_CGR_event.csv")
                + " --sv " + args.sv + " --out " + os.path.join(complex_path, "summary.png")):
            history.append('draw')
            update_history(history_file, 'draw')

        if args.genome_version == 'hg19' and args.gender == 'female' and cmd(
                "conda run --no-capture-output -n main Rscript "+ os.path.join(config.SCRIPT_DIR,"draw","ccrr_draw.R")+" --ref "+ os.path.join(config.SCRIPT_DIR,"draw/hg19_wm.txt")+" --cn "
                + args.cn + " --shatterseek " + os.path.join(complex_path, "shatterseek", "chromothripsis_summary.csv")
                + " --region " + os.path.join(config.SCRIPT_DIR,"draw/region/hg19_wm.region")
                + " --ctlpscanner " + os.path.join(complex_path, "ctlpscanner", "CTLPRegion.txt")
                + " --sa " + os.path.join(complex_path, "SA", "SA_amplicons.csv")
                + " --aa " + os.path.join(complex_path, "AA", args.tumor_id + "_classification",args.tumor_id + "_result_data.json")
                + " --gGnome " + os.path.join(complex_path, "gGnome", "event_footprints.txt")
                + " --starfish " + os.path.join(complex_path, "starfish", args.tumor_id + "_connected_CGR_event.csv")
                + " --sv " + args.sv + " --out " + os.path.join(complex_path, "summary.png")):
            
            history.append('draw')
            update_history(history_file, 'draw')

        if args.genome_version == 'hg38' and args.gender == 'male' and cmd(
                "conda run --no-capture-output -n main Rscript "+ os.path.join(config.SCRIPT_DIR,"draw","ccrr_draw.R")+" --ref "+ os.path.join(config.SCRIPT_DIR,"draw/hg38_man.txt") + " --cn "
                + args.cn + " --shatterseek " + os.path.join(complex_path, "shatterseek", "chromothripsis_summary.csv")
                + " --region " + os.path.join(config.SCRIPT_DIR,"draw/region/hg38_man.region")
                + " --ctlpscanner " + os.path.join(complex_path, "ctlpscanner", "CTLPRegion.txt")
                + " --sa " + os.path.join(complex_path, "SA", "SA_amplicons.csv")
                + " --aa " + os.path.join(complex_path, "AA", args.tumor_id + "_classification",args.tumor_id + "_result_data.json")
                + " --gGnome " + os.path.join(complex_path, "gGnome", "event_footprints.txt")
                + " --starfish " + os.path.join(complex_path, "starfish", args.tumor_id + "_connected_CGR_event.csv")
                + " --sv " + args.sv + " --out " + os.path.join(complex_path, "summary.png")):
            history.append('draw')
            update_history(history_file, 'draw')

        if args.genome_version == 'hg38' and args.gender == 'female' and cmd(
                "conda run --no-capture-output -n main Rscript "+ os.path.join(config.SCRIPT_DIR,"draw","ccrr_draw.R")+" --ref "+ os.path.join(config.SCRIPT_DIR,"draw/hg38_wm.txt")+" --cn "
                + args.cn + " --shatterseek " + os.path.join(complex_path, "shatterseek", "chromothripsis_summary.csv")
                + " --region " + os.path.join(config.SCRIPT_DIR,"draw/region/hg38_wm.region")
                + " --ctlpscanner " + os.path.join(complex_path, "ctlpscanner", "CTLPRegion.txt")
                + " --sa " + os.path.join(complex_path, "SA", "SA_amplicons.csv")
                + " --aa " + os.path.join(complex_path, "AA", args.tumor_id + "_classification",args.tumor_id + "_result_data.json")
                + " --gGnome " + os.path.join(complex_path, "gGnome", "event_footprints.txt")
                + " --starfish " + os.path.join(complex_path, "starfish", args.tumor_id + "_connected_CGR_event.csv")
                + " --sv " + args.sv + " --out " + os.path.join(complex_path, "summary.png")):
            history.append('draw')
            update_history(history_file, 'draw')


if __name__ == "__main__":
    main()

