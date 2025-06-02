import argparse
import os
import json
import csv

def validate_file(path):
    if not os.path.isfile(path):
        raise argparse.ArgumentTypeError(f"Invalid file path: {path}")
    return path

def save_json(data, path):
    with open(path, "w") as f:
        json.dump(data, f, indent=2)


def parse_int(x):
    x = x.strip()
    if x == "NA" or x == "":
        return None
    try:
        return int(x)
    except ValueError:
        return None

def parse_float(x):
    x = x.strip()
    if x == "NA" or x == "":
        return None
    try:
        return float(x)
    except ValueError:
        return None

def parse_arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument("--genome-version",choices=["hg19", "hg38"], required=True)    
    parser.add_argument("--region", type=validate_file, required=True)
    parser.add_argument("--centromere", type=validate_file, required=True)
    parser.add_argument("--cn", type=validate_file, required=True)
    parser.add_argument("--shatterseek", type=validate_file)
    parser.add_argument("--ctlpscanner", type=validate_file)
    parser.add_argument("--sa", type=validate_file)
    parser.add_argument("--aa", type=validate_file)
    parser.add_argument("--gGnome", type=validate_file)
    parser.add_argument("--starfish-events", type=validate_file)
    parser.add_argument("--starfish-class", type=validate_file)
    parser.add_argument("--sv", type=validate_file, required=True)
    parser.add_argument("--annotations", type=validate_file, required=True,)
    parser.add_argument("--output")

    return parser.parse_args()

def load_region(region_path):
    region_data = {}

    if not region_path or not os.path.exists(region_path) or os.path.getsize(region_path) == 0:
        return region_data

    with open(region_path) as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 3:
                continue
            try:
                chrom = fields[0]
                start = int(fields[1])
                end = int(fields[2])
                region_data[chrom] = {
                    "start": start,
                    "end": end
                }
            except (ValueError, IndexError):
                continue

    return region_data

def load_centromere(centromere_path):
    centromere_data = {}

    if not centromere_path or not os.path.exists(centromere_path) or os.path.getsize(centromere_path) == 0:
        return centromere_data

    with open(centromere_path) as f:
        for line in f:
            fields = line.strip().split()
            if len(fields) < 3:
                continue
            try:
                chrom = fields[0]
                length = int(fields[1])
                centromere_pos = int(fields[2])
                centromere_data[chrom] = {
                    "length": length,
                    "centromere": centromere_pos
                }
            except (ValueError, IndexError):
                continue

    return centromere_data

def load_cn(cn_path):
    cn_data = {}

    if not cn_path or not os.path.exists(cn_path) or os.path.getsize(cn_path) == 0:
        return cn_data

    with open(cn_path) as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 5:
                continue
            try:
                chrom = fields[0]
                if chrom == "chrY":
                    continue
                start = int(fields[1])
                end = int(fields[2])
                cn_int = int(fields[3])
                cn_float = float(fields[4])

                if chrom not in cn_data:
                    cn_data[chrom] = []

                cn_data[chrom].append({
                    "start": start,
                    "end": end,
                    "cn": cn_int,
                    "cn_float": cn_float
                })
            except (ValueError, IndexError):
                continue

    return cn_data

def load_sv(sv_path):
    sv_data = []
    if not sv_path or not os.path.exists(sv_path) or os.path.getsize(sv_path) == 0:
        return sv_data

    with open(sv_path) as f:
        try:
            header = next(f).rstrip("\n").split("\t")
        except StopIteration:
            return sv_data

        col_idx = {name: i for i, name in enumerate(header)}
        idx = col_idx.get
        for line in f:
            fields = line.rstrip("\n").split("\t")
            try:
                chrom1 = fields[idx("chrom1")] if idx("chrom1") is not None else None
                pos1   = fields[idx("pos1")]   if idx("pos1")   is not None else None
                chrom2 = fields[idx("chrom2")] if idx("chrom2") is not None else None
                pos2   = fields[idx("pos2")]   if idx("pos2")   is not None else None
                svtype = fields[idx("svtype")] if idx("svtype") is not None else None

            except IndexError:
                continue

            if None in (chrom1, pos1, chrom2, pos2, svtype):
                continue
            if chrom1 == "chrY" or chrom2 == "chrY":
                continue

            try:
                pos1 = int(pos1)
                pos2 = int(pos2)
            except ValueError:
                continue

            strand1 = (
                fields[idx("strand1")] if idx("strand1") is not None and idx("strand1") < len(fields) else "."
            )
            strand2 = (
                fields[idx("strand2")] if idx("strand2") is not None and idx("strand2") < len(fields) else "."
            )
            try:
                level = (
                    int(fields[idx("level")]) if idx("level") is not None and idx("level") < len(fields) else 0
                )
            except ValueError:
                level = 0
            try:
                evidence = (
                    fields[idx("from")] if idx("from") is not None and idx("from") < len(fields) else "."
                )
                if isinstance(evidence, str) and evidence.startswith("[") and evidence.endswith("]"):
                    evidence = evidence.strip("[]").replace("'", "").replace('"', '')
                    evidence = "; ".join([x.strip() for x in evidence.split(",") if x.strip()])
            except Exception:
                evidence = "."
            sv_data.append(
                {
                    "chrom1": chrom1,
                    "pos1": pos1,
                    "chrom2": chrom2,
                    "pos2": pos2,
                    "strand1": strand1,
                    "strand2": strand2,
                    "evidence":evidence,
                    "level": level,
                    "svtype": svtype,
                }
            )
    return sv_data

def load_shatterseek(shatterseek_path):
    data = []
    if not shatterseek_path or not os.path.exists(shatterseek_path) or os.path.getsize(shatterseek_path) == 0:
        return data
    with open(shatterseek_path, "r") as f:
        header = next(f)
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 1:
                continue

            target_cols = 27
            if len(fields) < target_cols:
                fields += ["NA"] * (target_cols - len(fields))
            else:
                fields = fields[:target_cols]

            chrom = fields[0].strip()
            start = parse_int(fields[1])
            end = parse_int(fields[2])
            number_DEL = parse_int(fields[3])
            number_DUP = parse_int(fields[4])
            number_h2hINV = parse_int(fields[5])
            number_t2tINV = parse_int(fields[6])
            number_TRA = parse_int(fields[7])
            clusterSize_including_TRA = parse_int(fields[8])
            number_SVs_sample = parse_int(fields[9])
            number_CNV_segments = parse_int(fields[10])
            pval_fragment_joins = parse_float(fields[11])
            chr_breakpoint_enrichment = parse_float(fields[12])
            pval_exp_chr = parse_float(fields[13])
            pval_exp_cluster = parse_float(fields[14])
            max_number_osc_2_states = parse_int(fields[15])
            max_number_osc_3_states = parse_int(fields[16])
            number_CN_segments_chr = parse_int(fields[17])
            max_number_osc_2_states_chr = parse_int(fields[18])
            max_number_osc_3_states_chr = parse_int(fields[19])
            inter_number_DEL = parse_int(fields[20])
            inter_number_h2hINV = parse_int(fields[21])
            inter_number_t2tINV = parse_int(fields[22])
            inter_number_DUP = parse_int(fields[23])
            inter_pval_fragment_joins = parse_float(fields[24])

            inter_other_chroms = fields[25].strip()
            if inter_other_chroms in ("NA", ""):
                inter_other_chroms = None

            inter_other_chroms_coords_all = fields[26].strip()
            if inter_other_chroms_coords_all in ("NA", ""):
                inter_other_chroms_coords_all = None

            if chrom and not chrom.startswith("chr"):
                chrom = "chr" + chrom

            if start is None or end is None:
                continue

            def safe_sum(*vals):
                return sum(v for v in vals if v is not None)

            sum1 = safe_sum(
                number_DEL, number_DUP, number_h2hINV, number_t2tINV,
                number_TRA, clusterSize_including_TRA
            )
            sum2 = safe_sum(
                inter_number_DEL, inter_number_h2hINV, inter_number_t2tINV, inter_number_DUP
            )
            max_osc_CN_seg = max_number_osc_2_states if max_number_osc_2_states is not None else 0

            classification = None
            if (sum1 >= 6 and max_osc_CN_seg >= 7) or (sum1 >= 3 and sum2 >= 4 and max_osc_CN_seg >= 7):
                classification = "high"
            elif sum1 >= 6 and max_osc_CN_seg >= 4:
                classification = "low"

            if classification is None:
                continue

            row_data = {
                "chrom": chrom if chrom != "NA" else None,
                "start": start,
                "end": end,
                "number_DEL": number_DEL,
                "number_DUP": number_DUP,
                "number_h2hINV": number_h2hINV,
                "number_t2tINV": number_t2tINV,
                "number_TRA": number_TRA,
                "clusterSize_including_TRA": clusterSize_including_TRA,
                "number_SVs_sample": number_SVs_sample,
                "number_CNV_segments": number_CNV_segments,
                "pval_fragment_joins": pval_fragment_joins,
                "chr_breakpoint_enrichment": chr_breakpoint_enrichment,
                "pval_exp_chr": pval_exp_chr,
                "pval_exp_cluster": pval_exp_cluster,
                "max_number_oscillating_CN_segments_2_states": max_number_osc_2_states,
                "max_number_oscillating_CN_segments_3_states": max_number_osc_3_states,
                "number_CN_segments_chr": number_CN_segments_chr,
                "max_number_oscillating_CN_segments_2_states_chr": max_number_osc_2_states_chr,
                "max_number_oscillating_CN_segments_3_states_chr": max_number_osc_3_states_chr,
                "inter_number_DEL": inter_number_DEL,
                "inter_number_h2hINV": inter_number_h2hINV,
                "inter_number_t2tINV": inter_number_t2tINV,
                "inter_number_DUP": inter_number_DUP,
                "inter_pval_fragment_joins": inter_pval_fragment_joins,
                "inter_other_chroms": inter_other_chroms,
                "inter_other_chroms_coords_all": inter_other_chroms_coords_all,
                "classification": classification
            }
            data.append(row_data)

    return data

def load_ctlpscanner(path):
    result = []
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return result
    with open(path, "r") as f:
        next(f)
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 10:
                continue
            try:
                win_no        = int(float(fields[1]))
                win_size      = float(fields[2])
                chrom_int     = int(float(fields[3]))
                start         = int(float(fields[4]))
                end           = int(float(fields[5]))
                exp_switch_no = float(fields[6])
                switch_no     = int(float(fields[7]))
                lr_val        = float(fields[8])
                pvalue_val    = float(fields[9])
            except:
                continue

            if chrom_int == 23 or chrom_int < 1 or chrom_int > 22:
                continue
            chrom_str = f"chr{chrom_int}" if chrom_int < 22 else "chrX"
            if start is None or end is None:
                continue

            entry = {
                "WinNo": win_no,
                "WinSize": win_size,
                "Chrom": chrom_str,
                "Start": start,
                "End": end,
                "ExpSwitchNo": exp_switch_no,
                "SwitchNo": switch_no,
                "LR": lr_val,
                "Pvalue": pvalue_val
            }
            result.append(entry)

    return result


def load_gGnome(gGnome_path):
    event_data = {}
    if not gGnome_path or not os.path.exists(gGnome_path) or os.path.getsize(gGnome_path) == 0:
        return event_data
    with open(gGnome_path) as f:
        header = next(f)
        for line in f:
            fields = line.strip().split()
            if len(fields) < 4:
                continue
            chrom, start, end, description = fields[0], parse_int(fields[1]), parse_int(fields[2]), fields[3]
            if chrom and not chrom.startswith("chr"):
                chrom = "chr" + chrom
            if start is None or end is None:
                continue
            if description not in event_data:
                event_data[description] = []
            event_data[description].append({
                "chrom": chrom,
                "start": start,
                "end": end
            })
    return event_data



def load_starfish(starfish_events_path, starfish_class_path):
    events = {}
    def safe_open_csv(path):
        return (
            open(path, newline='') if path and os.path.exists(path) and os.path.getsize(path) > 0 else None
        )

    f1 = safe_open_csv(starfish_events_path)
    if f1:
        reader = csv.DictReader(f1)
        for row in reader:
            chrom = row.get("chr", "").strip()
            if chrom and not chrom.startswith("chr"):
                chrom = "chr" + chrom
            event = {
                "chr": chrom,
                "start": parse_int(row.get("start", "")),
                "end": parse_int(row.get("end", "")),
                "sample": row.get("sample", "").strip(),
                "CGR_status": row.get("CGR_status", "").strip(),
                "link_chromosome": row.get("link_chromosome", "").strip(),
                "cluster_id": row.get("cluster_id", "").strip()
            }
            events[event["cluster_id"]] = event

    f2 = safe_open_csv(starfish_class_path)
    if f2:
        reader = csv.DictReader(f2)
        for row in reader:
            cluster_id = row.get("cluster_id", "").strip()
            classification_info = {
                "Brk_dispersion_MAD_mean_total": parse_float(row.get("Brk_dispersion_MAD_mean_total", "")),
                "Loss_size_percentage": parse_float(row.get("Loss_size_percentage", "")),
                "Gain_size_percentage": parse_float(row.get("Gain_size_percentage", "")),
                "log_max_CN": parse_float(row.get("log_max_CN", "")),
                "max_telo_loss_percentage": parse_float(row.get("max_telo_loss_percentage", "")),
                "CGR_signature": row.get("CGR_signature", "").strip()
            }
            if cluster_id in events:
                events[cluster_id].update(classification_info)
            else:
                classification_info["cluster_id"] = cluster_id
                events[cluster_id] = classification_info

    return events

def load_sa(sa_amplicons_path):
    results = []
    if not sa_amplicons_path or not os.path.exists(sa_amplicons_path) or os.path.getsize(sa_amplicons_path) == 0:
        return results
    with open(sa_amplicons_path, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            seqnames = row.get("seqnames", "").strip()
            entry = {
                "seqnames": seqnames,
                "start": parse_int(row.get("start", "")),
                "end": parse_int(row.get("end", "")),
                "width": parse_int(row.get("width", "")),
                "strand": row.get("strand", "").strip(),
                "id": parse_int(row.get("id", "")),
                "nSegments": parse_int(row.get("nSegments", "")),
                "medianCN": parse_int(row.get("medianCN", "")),
                "maxCN": parse_int(row.get("maxCN", "")),
                "size_amplicon": parse_int(row.get("size_amplicon", "")),
                "nChrs_amplicon": parse_int(row.get("nChrs_amplicon", "")),
                "nRegions_amplicon": parse_int(row.get("nRegions_amplicon", "")),
                "medianCN_amplicon": parse_int(row.get("medianCN_amplicon", "")),
                "cnSpan_amplicon": parse_int(row.get("cnSpan_amplicon", "")),
                "cnStates_amplicon": parse_int(row.get("cnStates_amplicon", "")),
                "nSegments_amplicon": parse_int(row.get("nSegments_amplicon", "")),
                "nSVs_amplicon": parse_int(row.get("nSVs_amplicon", "")),
                "nSVsInternal_amplicon": parse_int(row.get("nSVsInternal_amplicon", ""))
            }
            results.append(entry)
    return results

def load_aa(aa_path):
    if not aa_path or not os.path.exists(aa_path) or os.path.getsize(aa_path) == 0:
        return {}
    with open(aa_path, "r") as f:
        data = json.load(f)
    return data

def load_annotations(annotations_path):
    annotations = []
    if not annotations_path or not os.path.exists(annotations_path) or os.path.getsize(annotations_path) == 0:
        return annotations
    with open(annotations_path, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            gene_id = row.get("GeneId", "").strip()
            gene_name = row.get("GeneName", "").strip()
            chrom = row.get("Chromosome", "").strip()
            if chrom and not chrom.startswith("chr"):
                chrom = "chr" + chrom
            strand = row.get("Strand", "").strip()
            gene_start = parse_int(row.get("GeneStart", ""))
            gene_end = parse_int(row.get("GeneEnd", ""))
            karyotype_band = row.get("KaryotypeBand", "").strip()
            synonyms_field = row.get("Synonyms", "").strip()
            if synonyms_field in ("NA", ""):
                synonyms = []
            else:
                synonyms = [s.strip() for s in synonyms_field.split(';') if s.strip()]
            annotation = {
                "GeneId": gene_id,
                "GeneName": gene_name,
                "Chromosome": chrom,
                "Strand": strand,
                "GeneStart": gene_start,
                "GeneEnd": gene_end,
                "KaryotypeBand": karyotype_band,
                "Synonyms": synonyms
            }
            annotations.append(annotation)
    return annotations

def main():
    args = parse_arguments()
    region_data = load_region(args.region)
    centromere_data = load_centromere(args.centromere)
    cn_data = load_cn(args.cn)
    sv_data = load_sv(args.sv)
    shatterseek_data = load_shatterseek(args.shatterseek)
    ctlp_data = load_ctlpscanner(args.ctlpscanner)
    gGnome_data = load_gGnome(args.gGnome)
    starfish_data = load_starfish(args.starfish_events, args.starfish_class)
    sa_data = load_sa(args.sa)
    aa_data = load_aa(args.aa)
    annotations_data = load_annotations(args.annotations)
    output_data = {
        "region": region_data,
        "centromere": centromere_data,
        "cn": cn_data,
        "sv": sv_data,
        "shatterseek": shatterseek_data,
        "ctlpscanner": ctlp_data,
        "gGnome": gGnome_data,
        "starfish": starfish_data,
        "sa": sa_data,
        "aa": aa_data,
        "annotations": annotations_data
    }

    save_json(output_data, args.output)

    

if __name__ == "__main__":
    main()