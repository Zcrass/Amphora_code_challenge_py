import os
import gzip
import pandas as pd
# pd.options.mode.chained_assignment = None  # default='warn'

def get_vcf_names(vcf_path):
    ''' 
    This function identify the start of the vcf table and 
    return the names of it.
    Modified from: https://www.biostars.org/p/416324/#9480044
    Paremeters: 
        vcf_path: the VCF file to be read. Could be gziped

    Retrurns:
        A list of names taken from the VCF file

    '''
    with open(vcf_path, 'rb') as test_f:
        gziped = test_f.read(2) == b'\x1f\x8b'
    if gziped == True:
        opening_fun = gzip.open(vcf_path, 'rt')
    else:
        opening_fun = open(vcf_path, 'rt')
    with opening_fun as ifile:
        for line in ifile:
            if line.startswith('#CHROM'):
                vcf_names = [x for x in line.split('\t')]
                vcf_names[0] = vcf_names[0].replace('#', '')
                vcf_names[-1] = vcf_names[-1].rstrip()
                break
    ifile.close()
    return vcf_names

def import_geno(file, format):
    ''' 
    This function reads a CSV or VCF file from one individual sample
    and import the content as a data frame.
    VCF files should include ten columns: ['CHROM', 'POS', 'ID', 'REF',
    'ALT', 'QUAL','FILTER', 'INFO', 'FORMAT', 'individual genotye'].
    VCF files could be plain text or gziped.
    CSV files should include five columns: ['CHROM;POS', 'REF',
    'ALT', 'ALT.1', 'individual genotye'].
    
    Parameters:
        file: input file to be processed. Could be either a CSV or a VCF.

        format: either 'csv' or 'vcf'

    Returns:
        A dataframe of 5 columns: 
        ['CHROM', 'POS', 'REF', 'ALT', individual genotype].

    '''

    ind = os.path.basename(file).split('.', 1)[0]
    if format == 'vcf':
        names = get_vcf_names(file)
        with open(file, 'rb') as test_f:
            if test_f.read(2) == b'\x1f\x8b':
                compression = 'gzip'
            else: 
                compression = None
        geno = pd.read_csv(file, comment='#', compression=compression,
                             delim_whitespace=True, header=None,
                             names=names)
        geno = geno[['CHROM', 'POS', 'REF', 'ALT', ind]]
    elif format == 'csv':
        geno = pd.read_csv(file)
        geno = pd.concat([geno.iloc[:, 0].str.split(pat=';', expand=True),
                          geno[['REF', 'ALT', ind]]], axis=1) 
        geno.columns = ['CHROM', 'POS', 'REF', 'ALT', ind]
    geno['CHROM'] = geno['CHROM'].astype('int64')
    geno['POS'] = geno['POS'].astype('int64')
    geno[ind] = geno[ind].str.replace(pat=',', repl='|', regex=False)
    return geno

def geno2num(genotye):
    '''
    This function 'turn a series of genotype format values (0|0)
    into a series of numeric genotypes (0, 1, 2)
    Values of 0|0 are homozygous for the reference allele and are converted to 0
    Values of 0|1 or 1|0 are heterozygous and are converted to 1
    Values of 1|1 are homozygous for the alternative allele and are converted to 2
  
    Parameters:
        genotype: Series of genotype values in format '0|0', '0|1',  '1|1'
  
    Returns:
        A series of numeric values representing the genotypes (0, 1, 2)

    '''
    new_geno = genotye.str.replace('1|1', '2', regex=False)
    new_geno = new_geno.str.replace('1|0', '1', regex=False)
    new_geno = new_geno.str.replace('0|1', '1', regex=False)
    new_geno = new_geno.str.replace('0|0', '0', regex=False)
    return new_geno
