import pandas as pd
import argparse

def main():
    parser = argparse.ArgumentParser(description='Convert CN file to AA format')
    parser.add_argument('-input', type=str, help='Path to the input CNV file')
    parser.add_argument('-output', type=str, help='Path to the output file')

    args = parser.parse_args()

    df = pd.read_csv(args.input, sep='\t', header=None)
    df = df.iloc[:, :-1]
    df.columns = ['chr', 'start', 'end', 'copy_number']
    df.to_csv(args.output, sep='\t', index=False,header=False)

if __name__ == "__main__":
    main()

