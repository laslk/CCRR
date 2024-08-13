# chrom1 pos1 chrom2 pos2 strand1 strand2 svtype quality filter from id
# manta delly svaba gridss lumpy soreca
import argparse
import multiprocessing
import os
import re
import vcf
import pandas as pd
from collections import Counter


lumpy_default_filter = 1
chr_accepted = ["chr" + str(i) for i in range(1, 23)] + ["chrX", "chrY"]
strand_map = {
    "5to5": {
        "same_chrom": ["-", "-", "t2tINV"],
        "diff_chrom": ["-", "-", "TRA"]
    },
    "3to3": {
        "same_chrom": ["+", "+", "h2hINV"],
        "diff_chrom": ["+", "+", "TRA"]
    },
    "3to5": {
        "same_chrom": ["+", "-", "DEL"],
        "diff_chrom": ["+", "-", "TRA"]
    },
    "5to3": {
        "same_chrom": ["-", "+", "DUP"],
        "diff_chrom": ["-", "+", "TRA"]
    }
}

def worker_process(start, end, svall, threshold, process_num):
    id_relation_sub = []
    total_iterations = (end - start) * len(svall)
    current_iteration = 0
    last_percentage = 0
    for ind1 in range(start, end):
        line1 = svall.iloc[ind1]
        for _, line2 in svall.iterrows():
            current_iteration += 1

            current_percentage = int((current_iteration / total_iterations) * 100)
            if current_percentage >= last_percentage + 5:
                print(f"Process-{process_num} Progress: {current_percentage}%")
                last_percentage = current_percentage

            if sameline(line1, line2, threshold):
                id_relation_sub.append([line1['id'], line2['id']])
    return id_relation_sub
def parallel_id_relation(svall, threshold, num_processes):
    pool = multiprocessing.Pool(processes=num_processes)
    chunk_size = len(svall) // num_processes

    tasks = []
    for i in range(num_processes):
        start_index = i * chunk_size
        end_index = start_index + chunk_size if i != num_processes - 1 else len(svall)
        tasks.append((start_index, end_index, svall, threshold, i))

    results = pool.starmap(worker_process, tasks)

    id_relation = []
    for sublist in results:
        id_relation.extend(sublist)

    pool.close()
    pool.join()

    return id_relation


class UnionFind:
    def __init__(self):
        self.parent = {}

    def find(self, node):
        if node not in self.parent:
            self.parent[node] = node
        elif self.parent[node] != node:
            self.parent[node] = self.find(self.parent[node])
        return self.parent[node]

    def union(self, node1, node2):
        root1 = self.find(node1)
        root2 = self.find(node2)

        if root1 != root2:
            self.parent[root1] = root2


def group_related_elements(relations):
    uf = UnionFind()

    for a, b in relations:
        uf.union(a, b)

    groups = {}
    for a, b in relations:
        root = uf.find(a)
        if root not in groups:
            groups[root] = set()
        groups[root].add(a)
        groups[root].add(b)

    return [list(group) for group in groups.values()]


def mergeline(id, id_svall):
    if type(id) == str:
        line = id_svall[id]
        new = [line['chrom1'], line['pos1'], line['chrom2'], line['pos2'], line['strand1'], line['strand2'],
               line['from'], line['id'], 1]
    elif type(id) == list:
        lines = [id_svall[i].tolist() for i in id]
        chrom1 = merge_chr_str([line[0] for line in lines])
        chrom2 = merge_chr_str([line[2] for line in lines])
        pos1 = merge_pos([line[1] for line in lines])
        pos2 = merge_pos([line[3] for line in lines])
        strand1 = merge_chr_str([line[4] for line in lines])
        strand2 = merge_chr_str([line[5] for line in lines])
        froms = list(set(line[9] for line in lines))
        ids = [line[10] for line in lines]
        new = [chrom1, pos1, chrom2, pos2, strand1, strand2, froms, ids, len(froms)]
    else:
        new = False
    return new


def merge_chr_str(items):
    items = [item for item in items if item != "None"]
    counts = Counter(items).most_common()
    if counts:
        return counts[0][0]
    else:
        return 'None'

def merge_pos(poses):
    poses = [item for item in poses if item != 'None']
    counts = Counter(poses).most_common()
    modes = [int(item) for item, count in counts if count == counts[0][1]]
    return round(sum(modes) / len(modes))


def chromosome_key(item):
    chrom_mapping = {'chrX': 'chr23', 'chrY': 'chr24'}
    chrom = item[0] if item[0] not in chrom_mapping else chrom_mapping[item[0]]
    pos = int(item[1])
    return chrom, pos


def validate_file(filepath):
    if not os.path.exists(filepath):
        raise argparse.ArgumentTypeError(f"The file '{filepath}' does not exist!")
    return filepath


def getarg():
    parser = argparse.ArgumentParser(description="merge SVs")

    parser.add_argument('-manta', type=validate_file, help="manta vcf result")
    parser.add_argument('-delly', type=validate_file, help="delly vcf result")
    parser.add_argument('-svaba', type=validate_file, help="svaba vcf result")
    parser.add_argument('-gridss', type=validate_file, help="GRIDSS vcf result")
    parser.add_argument('-lumpy', type=validate_file, help="LUMPY vcf result")
    parser.add_argument('-soreca', type=validate_file, help="soreca result")

    parser.add_argument('--manta-filter', type=int, default=None, help="Filter threshold for manta")
    parser.add_argument('--delly-filter', type=int, default=None, help="Filter threshold for delly")
    parser.add_argument('--svaba-filter', type=int, default=None, help="Filter threshold for svaba")
    parser.add_argument('--gridss-filter', type=int, default=None, help="Filter threshold for GRIDSS")
    parser.add_argument('--lumpy-filter', type=int, default=None, help="Filter threshold for LUMPY")

    parser.add_argument('--threshold', type=int, default=150, help="threshold for determination, defaults to 150bp")
    parser.add_argument('-t', type=int, default=8, help="Set the number of processes")
    merge_choices = ["intersection", "union", "x-or-more"]
    parser.add_argument('--merge-method', type=str, required=True, choices=merge_choices,
                        help="Choose a merging method: \n"
                             "1. 'intersection': Merges only the SVs that are identified by all SV callers.\n"
                             "2. 'union': Merges all SVs identified by any of the SV callers.\n"
                             "3. 'x-or-more': Merges SVs that are identified by at least x SV callers.\n"
                             "if only one svcaller is prepared , then this parameter is irrelevant")
    parser.add_argument('--primary-caller', type=str, choices=['None','manta', 'delly', 'svaba', 'gridss', 'lumpy', 'soreca'],default=None,
                        help="Specify the primary SV caller to keep all of its result.")

    parser.add_argument('-x', type=int, choices=[1, 2, 3, 4, 5, 6],
                        help="Specify the x. This argument is required when '--merge-method' is set to "
                             "'x-or-more'. Must be among the provided input files.")
    parser.add_argument('-o', type=str, default=".", help="output path")

    args = parser.parse_args()
    if args.merge_method == "x-or-more" and args.x is None:
        parser.error("--x is required when '--merge-method' is set to 'x-or-more'")
        raise SystemExit

    input = {
        'manta': args.manta, 'delly': args.delly, 'svaba': args.svaba, 'gridss': args.gridss, 'lumpy': args.lumpy,
        'soreca': args.soreca}
    filters = {
        'manta': args.manta_filter, 'delly': args.delly_filter, 'svaba': args.svaba_filter,
        'gridss': args.gridss_filter, 'lumpy': args.lumpy_filter}
    if args.primary_caller != None and args.primary_caller != 'None' and args.primary_caller not in input:
        parser.error(f"You must provide an input file for the primary SV caller '{args.primary_caller}'")
        raise SystemExit

    input_svfile = {k: v for k, v in input.items() if v is not None}
    merge_method = [args.merge_method, args.primary_caller, args.x]
    return input_svfile, merge_method, args.threshold, filters, args.o,args.t


def sameline(line1, line2, threshold):
    line1_chr1 = line1['chrom1'] if line1['chrom1'] != "None" else None
    line1_chr2 = line1['chrom2'] if line1['chrom2'] != "None" else None
    line1_pos1 = int(line1['pos1']) if line1['pos1'] != "None" else None
    line1_pos2 = int(line1['pos2']) if line1['pos2'] != "None" else None
    line1_str1 = line1['strand1'] if line1['strand1'] != "None" else None
    line1_str2 = line1['strand2'] if line1['strand2'] != "None" else None

    line2_chr1 = line2['chrom1'] if line2['chrom1'] != "None" else None
    line2_chr2 = line2['chrom2'] if line2['chrom2'] != "None" else None
    line2_pos1 = int(line2['pos1']) if line2['pos1'] != "None" else None
    line2_pos2 = int(line2['pos2']) if line2['pos2'] != "None" else None
    line2_str1 = line2['strand1'] if line2['strand1'] != "None" else None
    line2_str2 = line2['strand2'] if line2['strand2'] != "None" else None
    if line1['from'] == line2['from']:
        return False
    if line1_chr1 and not line1_chr2 and line2_chr1:
        if line1_chr1 == line2_chr1 and abs(line1_pos1 - line2_pos1) < threshold:
            if None in [line1_str1, line2_str1]:
                return True
            elif line1_str1 == line2_str1:
                return True
    elif not line1_chr1 and line1_chr2 and line2_chr2:
        if line1_chr2 == line2_chr2 and abs(line1_pos2 - line2_pos2) < threshold:
            if None in [line1_str2, line2_str2]:
                return True
            elif line1_str2 == line2_str2:
                return True
    elif line1_chr1 and line2_chr1 and line1_chr2 and line2_chr2:
        if line1_chr1 == line2_chr1 and line1_chr2 == line2_chr2:
            if abs(line1_pos1 - line2_pos1) < threshold and abs(line1_pos2 - line2_pos2) < threshold:
                if None in [line1_str1, line2_str1, line1_str2, line2_str2]:
                    return True
                elif line1_str1 == line2_str1 and line1_str2 == line2_str2:
                    return True
    return False


def get_strand(record1, record2):
    patterns = {
        "1": r"[a-zA-Z]+?\[\w+:\w+\[",
        "2": r"[a-zA-Z]+?\]\w+:\w+\]",
        "3": r"\[\w+:\w+\[[a-zA-Z]+",
        "4": r"\]\w+:\w+\][a-zA-Z]+"
    }
    alt1 = str(record1.ALT[0])
    alt2 = str(record2.ALT[0])
    if re.fullmatch(patterns["1"], alt1) and re.fullmatch(patterns["4"], alt2):
        return "3to5"
    elif re.fullmatch(patterns["4"], alt1) and re.fullmatch(patterns["1"], alt2):
        return "5to3"
    elif re.fullmatch(patterns["3"], alt1) and re.fullmatch(patterns["3"], alt2):
        return "5to5"
    elif re.fullmatch(patterns["2"], alt1) and re.fullmatch(patterns["2"], alt2):
        return "3to3"
    else:
        print(record1.ID + " strand error")
        return False


def delly_sv2bed(delly_vcf, filter):
    delly_sv = []

    strand_map = {
        "5to5": ["-", "-", "t2tINV"],
        "3to3": ["+", "+", "h2hINV"],
        "5to3": ["-", "+", "DUP"],
        "3to5": ["+", "-", "DEL"]
    }

    for record in vcf.Reader(open(delly_vcf, "r")):
        record_filter = "PASS" if not record.FILTER else record.FILTER[0]

        if filter is None or record.QUAL > filter:
            strand_info = strand_map.get(record.INFO['CT'])
            if not strand_info:
                continue

            if record.ID.startswith("BND") and record.CHROM in chr_accepted and record.INFO["CHR2"] in chr_accepted:
                delly_sv.append(
                    [record.CHROM, str(record.POS), record.INFO["CHR2"], str(record.INFO["POS2"]),
                     strand_info[0], strand_info[1], "TRA", record.QUAL, record_filter, "delly", "delly:" + record.ID]
                )

            elif record.CHROM in chr_accepted:
                delly_sv.append(
                    [record.CHROM, str(record.POS), record.CHROM, str(record.INFO['END']),
                     strand_info[0], strand_info[1], strand_info[2], record.QUAL, record_filter, "delly",
                     "delly:" + record.ID]
                )
    print(f"get {len(delly_sv)}  SV from delly")
    return delly_sv


def manta_sv2bed(manta_vcf, filter):
    manta_sv = []
    manta_dic = {record.ID: record for record in vcf.Reader(open(manta_vcf, "r"))}

    history = set()
    for id1, record1 in manta_dic.items():
        if id1 in history:
            continue

        record_filter = "PASS" if not record1.FILTER else record1.FILTER[0]
        if filter is None or record1.QUAL > filter:

            if record1.INFO["SVTYPE"] == "BND":
                id2 = str(record1.INFO['MATEID'][0])
                history.update([id1, id2])
                if id2 not in manta_dic:
                    continue
                record2 = manta_dic[id2]
                strand = get_strand(record1, record2)
                if strand is False:
                    continue
                chrom_relation = "same_chrom" if record1.CHROM == record2.CHROM else "diff_chrom"
                if record1.CHROM in chr_accepted and (
                        chrom_relation == "diff_chrom" or record2.CHROM in chr_accepted):
                    strand_info = strand_map[strand][chrom_relation]
                    manta_sv.append([record1.CHROM, str(record1.POS), record2.CHROM, str(record2.POS)] + strand_info
                                    + [record1.QUAL, record_filter, "manta", id1 + "-" + id2])

            else:
                history.add(id1)
                if record1.INFO["SVTYPE"] == "DEL":
                    manta_sv.append([record1.CHROM, str(record1.POS), record1.CHROM, str(record1.INFO["END"]), "+", "-",
                                     "DEL", record1.QUAL, record_filter, "manta", id1])
                elif record1.INFO["SVTYPE"] == "DUP":
                    manta_sv.append([record1.CHROM, str(record1.POS), record1.CHROM, str(record1.INFO["END"]), "-", "+",
                                     "DUP", record1.QUAL, record_filter, "manta", id1])
                else:
                    manta_sv.append([record1.CHROM, str(record1.POS), record1.CHROM, str(record1.INFO["END"]), "+", "-",
                                     "INS", record1.QUAL, record_filter, "manta", id1])

    print(f"get {len(manta_sv)}  from manta")
    return manta_sv


def lumpy_sv2bed(lumpy_vcf, filter):
    lumpy_sv = []
    lumpy_dic = {record.ID: record for record in vcf.Reader(open(lumpy_vcf, "r"))}

    if filter is None:
        sorted_ids = sorted(lumpy_dic.keys(), key=lambda x: lumpy_dic[x].QUAL, reverse=True)
        lumpy_filted = sorted_ids[:int(len(sorted_ids) * lumpy_default_filter)]
    else:
        lumpy_filted = [id for id, record in lumpy_dic.items() if record.QUAL > filter]

    print(f"get {len(lumpy_filted)} of {len(lumpy_dic)} SV from lumpy")

    history = set()

    for id1 in lumpy_filted:
        record1 = lumpy_dic[id1]
        if id1 in history or record1.CHROM not in chr_accepted:
            continue

        record_filter = "PASS" if not record1.FILTER else record1.FILTER[0]

        if record1.INFO["SVTYPE"] == "BND":
            id2 = record1.INFO["MATEID"][0]
            history.update([id1, id2])

            record2 = lumpy_dic.get(id2, None)
            if not record2:
                continue

            strand = get_strand(record1, record2)
            if strand is False:
                continue

            chrom_relation = "same_chrom" if record1.CHROM == record2.CHROM else "diff_chrom"
            if record1.CHROM in chr_accepted and (chrom_relation == "diff_chrom" or record2.CHROM in chr_accepted):
                strand_info = strand_map[strand][chrom_relation]
                lumpy_sv.append([record1.CHROM, str(record1.POS), record2.CHROM, str(record2.POS)] + strand_info
                                + [record1.QUAL, record_filter, "lumpy", "lumpy:" + id1 + "-" + id2])

        else:
            sv_type = record1.INFO["SVTYPE"]
            if sv_type == "INV":
                strand_a = record1.INFO["STRANDS"][0].split(":")[1]
                strand_b = record1.INFO["STRANDS"][1].split(":")[1]
                strand = "3to3" if strand_a >= strand_b else "5to5"
            else:
                strand = "3to5" if sv_type == "DEL" else "5to3"

            strand_info = strand_map[strand]["same_chrom"]
            lumpy_sv.append([record1.CHROM, str(record1.POS), record1.CHROM, str(record1.INFO["END"])] + strand_info
                            + [record1.QUAL, record_filter, "lumpy", "lumpy:" + record1.ID])

    return lumpy_sv


def svaba_sv2bed(svaba_vcf, filter):
    svaba_sv = []
    svaba_dic = {record.ID: record for record in vcf.Reader(open(svaba_vcf, "r"))}
    history = set()

    for id1, record1 in svaba_dic.items():
        record1 = svaba_dic[id1]
        record_filter = "PASS" if record1.FILTER == [] else record1.FILTER[0]
        id2 = record1.INFO["MATEID"]
        if id2 not in svaba_dic:
            continue
        record2 = svaba_dic[id2]
        if (filter is None or record1.QUAL > filter) and id1 not in history and id2 not in history:
            history.update([id1, id2])
            strand = get_strand(record1, record2)
            if strand is False:
                continue
            chrom_relation = "same_chrom" if record1.CHROM == record2.CHROM else "diff_chrom"
            strand_info = strand_map[strand][chrom_relation]
            svaba_sv.append(
                [record1.CHROM, str(record1.POS), record2.CHROM, str(record2.POS),
                 strand_info[0], strand_info[1], strand_info[2],
                 record1.QUAL, record_filter, "svaba", "svaba:" + id1 + "-" + id2]
            )
    print(f"get {len(svaba_sv)}  from svaba")
    return svaba_sv


def gridss_sv2bed(gridss_vcf, filter):
    gridss_sv = []
    gridss_dic = {record.ID: record for record in vcf.Reader(open(gridss_vcf, "r"))}

    history = set()

    for id1, record1 in gridss_dic.items():
        if id1 in history:
            continue

        record_filter = "PASS" if not record1.FILTER else record1.FILTER[0]
        if filter is not None and record1.QUAL <= filter:
            continue

        if record1.INFO["EVENTTYPE"] != "SGL":
            id2 = record1.INFO["MATEID"][0]
            record2 = gridss_dic[id2]
            history.update([id1, id2])

            strand = get_strand(record1, record2)
            if strand is False:
                continue

            chrom_relation = "same_chrom" if record1.CHROM == record2.CHROM else "diff_chrom"
            if record1.CHROM in chr_accepted and record2.CHROM in chr_accepted:
                strand_info = strand_map[strand][chrom_relation]
                gridss_sv.append([record1.CHROM, str(record1.POS), record2.CHROM, str(record2.POS)]
                                 + strand_info + [record1.QUAL, record_filter, "gridss", id1 + "-" + id2])
        else:
            history.add(id1)
            if record1.CHROM not in chr_accepted:
                continue
            alt_str = str(record1.ALT[0])
            if alt_str.startswith("."):
                gridss_sv.append([record1.CHROM, str(record1.POS), "None", "None", "+", "None", "SGL",
                                  record1.QUAL, record_filter, "gridss", id1 + '_1'])
                gridss_sv.append(["None", "None", record1.CHROM, str(record1.POS), "None", "+", "SGL",
                                  record1.QUAL, record_filter, "gridss", id1 + '_2'])
            elif alt_str.endswith("."):
                gridss_sv.append([record1.CHROM, str(record1.POS), "None", "None", "-", "None", "SGL",
                                  record1.QUAL, record_filter, "gridss", id1 + '_1'])
                gridss_sv.append(["None", "None", record1.CHROM, str(record1.POS), "None", "-", "SGL",
                                  record1.QUAL, record_filter, "gridss", id1 + '_2'])
    print(f"get {len(gridss_sv)}  from gridss")
    return gridss_sv


def read_soreca(soreca_unsnarl):
    soreca = pd.read_csv(soreca_unsnarl, sep=r"\s+", engine='python')
    soreca_sv = []
    for index, col in soreca.iterrows():
        soreca_sv.append([col["Exact_Loc_Pair_1"].split(":")[0], col["Exact_Loc_Pair_1"].split(":")[1],
                          col["Exact_Loc_Pair_2"].split(":")[0], col["Exact_Loc_Pair_2"].split(":")[1], "None", "None",
                          "None", "None", "None", "soreca", "soreca" + str(index)])

    print(f"get {len(soreca_sv)}  from soreca")
    return soreca_sv

def determine_svtype(row):
    if row['chrom1'] != row['chrom2']:
        return 'TRA'
    else:
        if row['strand1'] == '+' and row['strand2'] == '-':
            return 'DEL'
        elif row['strand1'] == '-' and row['strand2'] == '+':
            return 'DUP'
        elif row['strand1'] == '+' and row['strand2'] == '+':
            return 'h2hINV'
        elif row['strand1'] == '-' and row['strand2'] == '-':
            return 't2tINV'


def svmerge():
    input_svfile, merge_method, threshold, filters, output_path,processes = getarg()
    columns = ["chrom1", "pos1", "chrom2", "pos2", "strand1", "strand2", "svtype", "quality",
               "filter", "from", "id"]
    tools_mapping = {
        "delly": (delly_sv2bed, True),
        "manta": (manta_sv2bed, True),
        "lumpy": (lumpy_sv2bed, True),
        "svaba": (svaba_sv2bed, True),
        "gridss": (gridss_sv2bed, True),
        "soreca": (read_soreca, False)
    }

    svall = []
    for tool, (converter, uses_filter) in tools_mapping.items():
        if tool in input_svfile:
            if uses_filter:
                sv_data = pd.DataFrame(converter(input_svfile[tool], filter=filters.get(tool, None)), columns=columns)
            else:
                sv_data = pd.DataFrame(converter(input_svfile[tool]), columns=columns)

            svall.append(sv_data)

    if svall:
        svall = pd.concat(svall, axis=0)
    else:
        svall = pd.DataFrame(columns=columns)
        
    if len(input_svfile) > 1:
        id_svall = {line['id']: line for _, line in svall.iterrows()}
        id_relation = parallel_id_relation(svall, threshold, processes)
        id_relation = group_related_elements(id_relation)
        id_single = list(set(id_svall.keys()) - set({id for group in id_relation for id in group}))
        id_merged = id_single + id_relation
        merged_lines = [mergeline(id, id_svall) for id in id_merged
                        if type(id) == str and id_svall[id]['chrom1'] != "None"] + \
                       [mergeline(id, id_svall) for id in id_merged
                        if type(id) == list]
        method_map = {
            'intersection': len(input_svfile),
            'union': 1,
            'x-or-more': merge_method[2],
        }
        washed_lines = [line for line in merged_lines if line[8] >= method_map[merge_method[0]]]
        if merge_method[1]:
            washed_lines_primary = [line for line in merged_lines if merge_method[1] in line[6]]
            washed_lines = washed_lines + washed_lines_primary
        result_lines = []
        for line in washed_lines:
            if line not in result_lines:
                result_lines.append(line)
        result_lines = sorted(result_lines, key=chromosome_key)
    else:
        result_lines = [
            [line['chrom1'], line['pos1'], line['chrom2'], line['pos2'], line['strand1'], line['strand2'], line['from'], line['id'], 1] for _, line in svall.iterrows() if line['chrom1'] in chr_accepted and line['chrom2'] in chr_accepted and line['pos1'] != 'None' and line['pos2'] != 'None']
    sv_merged = pd.DataFrame(result_lines,
                             columns=["chrom1", "pos1", "chrom2", "pos2", "strand1", "strand2", "from", "id",
                                      "level"])
    filtered_data = sv_merged[(sv_merged['chrom1'].isin(chr_accepted)) &
                              (sv_merged['chrom2'].isin(chr_accepted)) &
                              (sv_merged['strand1'] != 'None') &
                              (sv_merged['strand2'] != 'None') &
                              (sv_merged['pos1'] != 'None') &
                              (sv_merged['pos2'] != 'None')]
    filtered_data = filtered_data[filtered_data['pos1'].astype(int) < filtered_data['pos2'].astype(int)]
    filtered_data.loc[:, 'svtype'] = filtered_data.apply(determine_svtype, axis=1)
    filtered_data.to_csv(os.path.join(output_path, "sv_merged.bed"), sep='\t', header=True, index=False)


if __name__ == "__main__":
    svmerge()
    print("*******************")
    print("***sv merge done***")
    print("*******************")
