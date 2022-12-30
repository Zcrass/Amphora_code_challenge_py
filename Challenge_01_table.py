#!/usr/bin/env python
import argparse
import logging as lg
import os
import pandas as pd
import sys
from functions import get_vcf_names, import_geno, geno2num


if __name__ == "__main__":
    lg.basicConfig(filename='Challenge_01.log', filemode='w',
                   format='%(name)s - %(levelname)s - %(message)s')
    logger = lg.getLogger('Challenge_01')
    logger.setLevel(lg.INFO)
    
    stdout_handler = lg.StreamHandler(sys.stdout)
    stdout_handler.setLevel(lg.INFO)
    logger.addHandler(stdout_handler)
    
    parser = argparse.ArgumentParser(prog = 'Challenge_01',
                                     description = 'Program to perform the first task in the challenge we use python to merge multiple genotype files in different format (VCF and CSV) from the Challenge samples folder into one large genotype table',
                                     epilog = "NANANA")
    parser.add_argument('-s', '--samples_path')
    parser.add_argument('-o', '--out_file')
    args = parser.parse_args()
    
    if args.samples_path == None:
        logger.error('ERROR: samples path is required!!!')
    elif args.out_file == None:
        logger.error('ERROR: out file is required!!!')
    else:  
        logger.info(f'listing files in folder {args.samples_path} ...')
        file_list = os.listdir(args.samples_path)
        vcf_files = [s for s in file_list if "vcf" in s.lower()]
        csv_files = [s for s in file_list if "csv" in s.lower()]
        logger.info(f'proccessing {len(file_list)} files')

        geno_table = pd.DataFrame(columns=["CHROM", "POS", "REF", "ALT",])
        j = 0
        for i in vcf_files:
            sample_file = args.samples_path + i
            j += 1
            print(j)
            sample_geno = import_geno(sample_file, format="vcf")
            # logger.info(f'')
            sample_geno.iloc[:, -1] = geno2num(sample_geno.iloc[:, -1])
            if geno_table.empty:
                geno_table = sample_geno
            else:
                geno_table = geno_table.merge(sample_geno, how="left")
        for i in csv_files:
            sample_file = args.samples_path + i
            j += 1
            print(j)
            sample_geno = import_geno(sample_file, format="csv")
            sample_geno.iloc[:, -1] = geno2num(sample_geno.iloc[:, -1])
            if geno_table.empty:
                geno_table = sample_geno
            else:
                geno_table = geno_table.merge(sample_geno, how="left")
        geno_table.to_csv((args.out_file), index=False)
        # print(geno_table.shape)


