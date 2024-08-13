import pandas as pd
import numpy as np
import argparse

def calculate_log2_ratio(input_file_path, output_file_path, sample_name, reference_copy_number=2):
    cn_data = pd.read_csv(input_file_path, sep='\t', header=None)
    cn_data.iloc[:, 4] = cn_data.iloc[:, 4].replace(0, np.nan)
    cn_data['log2_ratio'] = np.log2(cn_data.iloc[:, 4] / reference_copy_number)
    cn_data['log2_ratio'] = cn_data['log2_ratio'].where(cn_data.iloc[:, 4].notna(), np.nan)
    cn_data.insert(0, 'sample_name', sample_name)

    cn_data.drop(cn_data.columns[4], axis=1, inplace=True)
    cn_data.drop(cn_data.columns[4], axis=1, inplace=True)
    cn_data.iloc[:, 1] = cn_data.iloc[:, 1].str.replace('chr', '')
    header = ['sample', 'chro', 'start', 'stop', 'mean']
    cn_data.to_csv(output_file_path, sep='\t', index=False, header=header)

def main():
    parser = argparse.ArgumentParser(description='Calculate log2 ratio for CNV data.')
    parser.add_argument('input_file', type=str, help='Path to the input CNV file')
    parser.add_argument('output_file', type=str, help='Path to the output file')
    parser.add_argument('--sample-name', type=str, default='Sample', help='Name of the sample')

    args = parser.parse_args()

    calculate_log2_ratio(args.input_file, args.output_file, args.sample_name)

if __name__ == "__main__":
    main()


