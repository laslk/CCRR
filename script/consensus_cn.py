import argparse
import os
import json
from decimal import Decimal, getcontext
from math import log

chr_accepted = ["chr" + str(i) for i in range(1, 23)] + ["chrX", "chrY"]
getcontext().prec = 10

def get_delly(path, limit):
    delly_bp_dict = {}
    delly_cn_dict = {}
    with open(path, "r") as delly_seg:
        for line in delly_seg:
            chr_name, start, end, _, cn = line.split("\t")

            if chr_name in chr_accepted:
                if chr_name not in delly_bp_dict:
                    delly_bp_dict[chr_name] = [1]
                    delly_cn_dict[chr_name] = {"1": 2}
                delly_bp_dict[chr_name].append(int(start))
                delly_bp_dict[chr_name].append(int(end))
                delly_cn_dict[chr_name][start] = cn.split("\n")[0]
                delly_cn_dict[chr_name][end] = 2

    for chr_name, bps in delly_bp_dict.items():
        bps.append(limit[chr_name])
    return delly_bp_dict, delly_cn_dict


def get_sclust(path):
    sclust_bp_dict = {}
    sclust_cn_dict = {}
    last_end = {}

    with open(path, "r") as sclust_seg:
        for line in sclust_seg:
            _, chr_name, start, end, _, cn = line.split("\t")

            if chr_name in chr_accepted:
                if chr_name not in sclust_bp_dict:
                    sclust_bp_dict[chr_name] = []
                    sclust_cn_dict[chr_name] = {}

                sclust_bp_dict[chr_name].append(int(start))
                sclust_cn_dict[chr_name][start] = cn.split("\n")[0]
                last_end[chr_name] = int(end)

        for chr_name, end_val in last_end.items():
            sclust_bp_dict[chr_name].append(end_val)

    return sclust_bp_dict, sclust_cn_dict


def get_purple(path):
    purple_bp_dict = {}
    purple_cn_dict = {}
    last_end = {}

    with open(path, "r") as purple_seg:
        for line in purple_seg:
            if line.startswith("chromosome"):
                continue
            chr_name, start, end, cn, *_ = line.split("\t")

            if chr_name in chr_accepted:
                if chr_name not in purple_bp_dict:
                    purple_bp_dict[chr_name] = []
                    purple_cn_dict[chr_name] = {}

                purple_bp_dict[chr_name].append(int(start))
                purple_cn_dict[chr_name][start] = cn.split("\n")[0]
                last_end[chr_name] = int(end)

        for chr_name, end in last_end.items():
            purple_bp_dict[chr_name].append(end)
    return purple_bp_dict, purple_cn_dict


def get_cnvkit(path):
    cnvkit_bp_dict = {}
    cnvkit_cn_dict = {}
    last_end = {}

    with open(path, "r") as cnvkit_seg:
        for line in cnvkit_seg:
            chr_name, start, end, _, cn = line.split("\t")

            if chr_name in chr_accepted:
                if chr_name not in cnvkit_bp_dict:
                    cnvkit_bp_dict[chr_name] = [1, ]
                    cnvkit_cn_dict[chr_name] = {"1": 2}

                cnvkit_bp_dict[chr_name].append(int(start))
                cnvkit_cn_dict[chr_name][start] = cn.split("\n")[0]
                last_end[chr_name] = int(end)

        for chr_name, end in last_end.items():
            cnvkit_bp_dict[chr_name].append(end)
    return cnvkit_bp_dict, cnvkit_cn_dict

def get_sequenza(path):
    sequenza_bp_dict = {}
    sequenza_cn_dict = {}
    last_end = {}

    with open(path, "r") as sequenza_seg:
        for line in sequenza_seg:
            if line.startswith("chromosome"):
                continue
            chr_name, start,end,_,_,_,_,_,_, cn, *_ = line.split("\t")

            if chr_name in chr_accepted:
                if chr_name not in sequenza_bp_dict:
                    sequenza_bp_dict[chr_name] = [1]
                    sequenza_cn_dict[chr_name] = {"1": 2}

                sequenza_bp_dict[chr_name].append(int(start))
                sequenza_cn_dict[chr_name][start] = cn
                last_end[chr_name] = int(end)

        for chr_name, end in last_end.items():
            sequenza_bp_dict[chr_name].append(end)
    return sequenza_bp_dict, sequenza_cn_dict







def cn_smooth(bp_dict, cn_dict):
    for method, chrom_bps in bp_dict.items():
        for chrom, bps in chrom_bps.items():
            cn_values = [round(float(cn_dict[method][chrom][str(bp)])) for bp in bps[:-1]]
            smoothed_bps = [bps[0]] + [bps[i] for i in range(1, len(cn_values)) if cn_values[i] != cn_values[i - 1]] + [
                bps[-1]]
            smoothed_cn = {bp: round(float(cn_dict[method][chrom][str(bp)])) for bp in smoothed_bps[:-1]}
            bp_dict[method][chrom] = smoothed_bps
            cn_dict[method][chrom] = smoothed_cn

    return bp_dict, cn_dict


def get_breakpoint_regions(bp_dict, t, limit):
    regions_dict = {}

    for method, method_bps in bp_dict.items():
        chrom_region = {}

        for chrom, bps in method_bps.items():
            if not bps:
                continue

            def calculate_region(bps, i):
                left_limit = round((bps[i - 1] + bps[i]) / 2) if i > 0 else 1
                right_limit = round((bps[i] + bps[i + 1]) / 2) if i < len(bps) - 1 else limit[chrom]

                return [max(bps[i] - t, left_limit), min(bps[i] + t, right_limit)]

            regions = [calculate_region(bps, i) for i in range(len(bps))]
            chrom_region[chrom] = regions

        regions_dict[method] = chrom_region

    return regions_dict


def generate_boundaries(regions_dict, chrom):
    boundaries = set()
    for method, chrom_dict in regions_dict.items():
        if chrom in chrom_dict:
            for region in chrom_dict[chrom]:
                boundaries.add(region[0])
                boundaries.add(region[1])
    return sorted(list(boundaries))


def count_regions(boundaries, regions_dict, chrom):
    counts = []
    for i in range(len(boundaries) - 1):
        start = boundaries[i]
        end = boundaries[i + 1]
        count = 0
        for method, chrom_dict in regions_dict.items():
            if chrom in chrom_dict:
                for region in chrom_dict[chrom]:
                    if region[0] <= start < region[1] or region[0] < end <= region[1]:
                        count = count + 1
        counts.append(count)
    return counts


def remove_and_record_regions(segment_index, boundaries, regions_dict, chrom):
    regions_to_remove = []
    start = boundaries[segment_index]
    end = boundaries[segment_index + 1]
    for method in regions_dict:
        if chrom in regions_dict[method]:
            for region in regions_dict[method][chrom]:
                if region[0] <= start < region[1] or region[0] < end <= region[1]:
                    regions_to_remove.append(region)
                    regions_dict[method][chrom].remove(region)
    return regions_to_remove


def get_overlap_regions(regions_dict, chrom, target_count):
    recorded_regions = []
    while True:
        boundaries = generate_boundaries(regions_dict, chrom)
        counts = count_regions(boundaries, regions_dict, chrom)
        if target_count in counts:
            segment_index = counts.index(target_count)
            recorded_regions.append(
                remove_and_record_regions(segment_index, boundaries, regions_dict, chrom))
        else:
            break
    return recorded_regions


def get_confidence_bp(regions, bp_dict, chrom):
    min_boundary = min([region[0] for region in regions])
    max_boundary = max([region[1] for region in regions])

    common_start = max([region[0] for region in regions])
    common_end = min([region[1] for region in regions])

    candidate_bps = []
    for chrom_bps in bp_dict.values():
        candidate_bps.extend([bp for bp in chrom_bps.get(chrom, []) if min_boundary <= bp <= max_boundary])

    confidence_bp = min(candidate_bps, key=lambda bp: abs(bp - (common_start + common_end) / 2))

    return confidence_bp


def get_different_cns(region, methods_cns):
    start = region[0]
    end = region[1]
    methods_regions = {}
    for method, bp_dict in methods_cns.items():
        split_regions = []
        current_start = start
        closest_bp = max([bp for bp in bp_dict if bp <= current_start], default=start)
        current_copy_number = bp_dict.get(closest_bp)

        for bp, copy_number in sorted(bp_dict.items()):
            if start < bp < end:
                split_regions.append({(current_start, bp): current_copy_number})
                current_start = bp
                current_copy_number = copy_number
        split_regions.append({(current_start, end): current_copy_number})
        methods_regions[method] = split_regions
    return methods_regions


def get_consensus_cn(chrom, chrom_consensus_bps, cn_dict):
    method_cns = {}
    for method, chrom_cn in cn_dict.items():
        if chrom in chrom_cn:
            method_cns[method] = chrom_cn[chrom]
    consensus_cns = {}
    if not len(chrom_consensus_bps) < 2:
        for i in range(len(chrom_consensus_bps) - 1):
            start, end = chrom_consensus_bps[i], chrom_consensus_bps[i + 1]
            different_cns = get_different_cns([start, end], method_cns)
            length_dict = {}
            total_length = 0
            for k, v in different_cns.items():
                for sub_dict in v:
                    for tuple_key, value in sub_dict.items():
                        interval_length = tuple_key[1] - tuple_key[0] + 1
                        if value in length_dict:
                            length_dict[value] += interval_length
                        else:
                            length_dict[value] = interval_length
                        total_length += interval_length
            proportion_dict = {}
            for key, value in length_dict.items():
                proportion_dict[key] = value / total_length
            consensus_cn=0
            for cn, weight in proportion_dict.items():
                if  cn==None or weight == None:
                    print("Warning! Something wrong during CN merging!")
                else: 
                    consensus_cn = consensus_cn+cn*weight
            bias = {}
            volatility= {}
            for k, v in different_cns.items():
                tool_bias = Decimal(0)
                tool_volatility  = Decimal(0)
                for sub_dict in v:
                    for tuple_key, value in sub_dict.items():
                        interval_length = tuple_key[1] - tuple_key[0] + 1
                        percent = Decimal(interval_length) / Decimal(total_length)
                        tool_volatility=abs(Decimal(value) - Decimal(consensus_cn) / (Decimal(log(consensus_cn + 2)) + Decimal(1)) * Decimal(percent))
                        tool_bias += (Decimal(value) - Decimal(consensus_cn) / (Decimal(log(consensus_cn + 2)) + Decimal(1)) * Decimal(percent))
                bias[k] = float(tool_bias*1000) 
                volatility[k] =  float(tool_volatility*1000)
            consensus_cns[(start, end)] = (round(consensus_cn),consensus_cn,json.dumps(bias),json.dumps(volatility))
    return consensus_cns



def validate_file(filepath):
    if not os.path.exists(filepath):
        raise argparse.ArgumentTypeError(f"The file '{filepath}' does not exist!")
    return filepath


def consensus_cn():
    parser = argparse.ArgumentParser(description="merge cn")
    parser.add_argument('-sclust', type=validate_file, help="sclust cn result")
    parser.add_argument('-delly', type=validate_file, help="delly cn result")
    parser.add_argument('-purple', type=validate_file, help="purple cn result")
    parser.add_argument('-cnvkit', type=validate_file, help="cnvkit cn result")
    parser.add_argument('-sequenza', type=validate_file, help="sequenza cn result")
    parser.add_argument('--threshold', type=int, default=5000,
                            help="threshold for determination, defaults to 5000bp")
    parser.add_argument('-o', type=str, default=".", help="output path")
    parser.add_argument('-ref', type=str, help="hg19 or hg38")
    args = parser.parse_args()
    cnv_keys = ['sclust', 'delly', 'purple', 'cnvkit','sequenza']
    cnv_dic = {key: getattr(args, key) for key in cnv_keys if getattr(args, key)}
    threshold = args.threshold
    outputpath = args.o
    ref = args.ref    
    if ref == "hg19":
        limit = {'chr1': 249250621, 'chr2': 243199373, 'chr3': 198022430,
                 'chr4': 191154276, 'chr5': 180915260, 'chr6': 171115067,
                 'chr7': 159138663, 'chr8': 146364022, 'chr9': 141213431,
                 'chr10': 135534747, 'chr11': 135006516, 'chr12': 133851895,
                 'chr13': 115169878, 'chr14': 107349540, 'chr15': 102531392,
                 'chr16': 90354753, 'chr17': 81195210, 'chr18': 78077248,
                 'chr19': 59128983, 'chr20': 63025520, 'chr21': 48129895,
                 'chr22': 51304566, 'chrX': 155270560, 'chrY': 59373566}
    else:
        limit = {'chr1': 248956422, 'chr2': 242193529, 'chr3': 198295559,
                 'chr4': 190214555, 'chr5': 181538259, 'chr6': 170805979,
                 'chr7': 159345973, 'chr8': 145138636, 'chr9': 138394717,
                 'chr10': 133797422, 'chr11': 135086622, 'chr12': 133275309,
                 'chr13': 114364328, 'chr14': 107043718, 'chr15': 101991189,
                 'chr16': 90338345, 'chr17': 83257441, 'chr18': 80373285,
                 'chr19': 58617616, 'chr20': 64444167, 'chr21': 46709983,
                 'chr22': 50818468, 'chrX': 156040895, 'chrY': 57227415}
    bp_dict = {}
    cn_dict = {}
    for method, path in cnv_dic.items():
        if method == "delly":
            bp_dict[method], cn_dict[method] = get_delly(path, limit)
        elif method == "sclust":
            bp_dict[method], cn_dict[method] = get_sclust(path)
        elif method == "purple":
            bp_dict[method], cn_dict[method] = get_purple(path)
        elif method == "cnvkit":
            bp_dict[method], cn_dict[method] = get_cnvkit(path)
        elif method == "sequenza":
            bp_dict[method], cn_dict[method] = get_sequenza(path)

    bp_dict, cn_dict = cn_smooth(bp_dict, cn_dict)
    regions_dict = get_breakpoint_regions(bp_dict, t=threshold, limit=limit)

    with open(os.path.join(outputpath, "consensus_cn.bed"), "w") as outfile,open(os.path.join(outputpath, "consensus_cn_statistic.bed"), "w") as sumfile:
        for chrom in chr_accepted:
            regions = []
            chrom_consensus_bps = []
            if len(cnv_dic) >= 2:
                for target_count in [5, 4, 3, 2]:
                    regions.extend(get_overlap_regions(regions_dict, chrom, target_count))
                for consensus_region in regions:
                    chrom_consensus_bps.append(get_confidence_bp(consensus_region, bp_dict, chrom))
                chrom_consensus_bps = sorted(chrom_consensus_bps)
            elif chrom in bp_dict[list(cnv_dic.keys())[0]]:
                chrom_consensus_bps = bp_dict[list(cnv_dic.keys())[0]][chrom]
            consensus_cns = get_consensus_cn(chrom, chrom_consensus_bps, cn_dict)
            
            if consensus_cns:
                for (start, end), (cn, cn_noround,bias,stable), in consensus_cns.items():
                    outfile.write(f"{chrom}\t{start}\t{end - 1}\t{cn}\t{cn_noround}\n")
                    sumfile.write(f"{chrom}\t{start}\t{end - 1}\t{cn}\t{cn_noround}\t{end - start}\t{bias}\t{stable}\n")

if __name__ == "__main__":
    consensus_cn()
    print("*******************")
    print("***cn merge done***")
    print("*******************")

