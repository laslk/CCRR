import warnings
import pandas as pd
import sys, getopt
warnings.simplefilter("ignore")
def main(argv):
    inputfile = ''
    outputfile = ''
    try:
        opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
    except getopt.GetoptError:
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg
    return inputfile,outputfile

def Write_File(LIST,filename):
    output = open(filename, "w")
    for i in range(len(LIST)):
        output.write(LIST[i] + '\n')

def Add_File(LIST,filename):
    LIST.to_csv(filename,mode='a',index=0,header=0,sep='\t')

def getfile_To_sclust(path):
    f = open(path, 'r')
    df = pd.read_csv(f,sep='\t',header=None,comment='#')
    print("read end")
    df.columns =['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','NORMAL','TUMOR']
    CHR = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16',
           'chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY']
    df = df[df['CHROM'].isin(CHR)]
    df_NOMAL= df['NORMAL'].str.split(":",expand=True)
    df_TUMOR = df['TUMOR'].str.split(":",expand=True)
    df_NOMAL.columns=['GT','GQ','DP','RD','AD','VAF','DP4']
    df_TUMOR.columns=['GT','GQ','DP','RD','AD','VAF','DP4']
    df_TUMOR['VAF'] = df_TUMOR['VAF'].str.split('%').str[0]
    df_NOMAL['VAF'] = df_NOMAL['VAF'].str.split('%').str[0]
    df = df[['CHROM','POS','ID','REF','ALT','QUAL','FILTER']]

    df_TUMOR['VAF'] = df_TUMOR['VAF'].astype(float)/100
    df_TUMOR['DP'] = df_TUMOR['DP'].astype(int)
    df_TUMOR['RD'] = df_TUMOR['RD'].astype(int)
    df_TUMOR['AD'] = df_TUMOR['AD'].astype(int)
    df_NOMAL['VAF'] = df_NOMAL['VAF'].astype(float)/100
    df_NOMAL['DP'] = df_NOMAL['DP'].astype(int)
    df_NOMAL['RD'] = df_NOMAL['RD'].astype(int)
    df_NOMAL['AD'] = df_NOMAL['AD'].astype(int)
    answer = pd.DataFrame()
    answer['DP_c'] = df_TUMOR['RD']+df_TUMOR['AD']
    answer['AF_c'] = df_TUMOR['AD']/answer['DP_c']
    answer['DP_n'] = df_NOMAL['RD']+df_NOMAL['AD']
    answer['AF_n'] = df_NOMAL['AD']/answer['DP_n']
    df = df.join(answer)

    message = df[(df['DP_c']>14) & (df['AF_c']>0.1) ]
    message[['DP_c', 'AF_c', 'DP_n', 'AF_n']] = message[['DP_c', 'AF_c', 'DP_n', 'AF_n']].astype(str)

    message['INFO'] = str('DP=') + message['DP_c'] + str(';AF=') + message['AF_c'] + str(';DP_N=') + message[
        'DP_n'] + str(';AF_N=') + message['AF_n']
    message = message[['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']]
    return message

def head():
    head = ['##fileformat=VCFv4.1',
            '##source=ChangefromVarScan2',
            '##INFO=<ID=DP,Number=1,Type=Integer, Description="Read Depth Tumor">',
            '##INFO=<ID=DP_N,Number=1,Type=Integer, Description="Read Depth Normal">',
            '##INFO=<ID=AF,Number=A,Type=Float, Description="Allelic Frequency Tumor">',
            '##INFO=<ID=AF_N,Number=A,Type=Float, Description="Allelic Frequency Normal">',
            '##INFO=<ID=FR,Number=1,Type=Float, Description="Forward-Reverse Score">',
            '##INFO=<ID=TG,Number=1,Type=String, Description="Target Name (Genome Partition)">',
            '##INFO=<ID=DB,Number=0,Type=Flag, Description="dbSNP Membership">',
            '##mfilterParameters= -af 0.2 -fr 0.2 -rc 5',
            '#CHROM\tPOS\tID\tREF\tALT\tQUA\tFILTER\tINFO']
    return head
inputfile,outputfile=main(sys.argv[1:])
message = getfile_To_sclust(inputfile)
print(message)
Write_File(head(),outputfile)
Add_File(message,outputfile)

