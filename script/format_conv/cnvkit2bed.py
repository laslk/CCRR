import pandas as pd
import argparse


def cns2bed(input_cns_path, output_bed_path):
    df = pd.read_csv(input_cns_path, sep="\t")
    df['log2_transformed'] = 2 ** (df['log2'] + 1)
    df_bed = df[['chromosome', 'start', 'end', 'cn', 'log2_transformed']]
    df_bed.to_csv(output_bed_path, sep="\t", index=False, header=False)


def main():
    parser = argparse.ArgumentParser(description="Convert CNVkit .cns files to .bed format with calculated CN values.")
    parser.add_argument("input_file", help="Input .cns file path")
    parser.add_argument("output_file", help="Output .bed file path")
    
    args = parser.parse_args()
    
    cns2bed(args.input_file, args.output_file)

if __name__ == "__main__":
    main()