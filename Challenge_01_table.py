#!/usr/bin/env python
import argparse
import logging as lg
import os
import pandas as pd
import sys
from functions import import_geno, geno2num


if __name__ == '__main__':
    lg.basicConfig(filename='Challenge_01.log', filemode='w', format='%(name)s - %(levelname)s - %(message)s')
    logger = lg.getLogger('Challenge_01')
    logger.setLevel(lg.INFO)
    
    stdout_handler = lg.StreamHandler(sys.stdout)
    stdout_handler.setLevel(lg.INFO)
    logger.addHandler(stdout_handler)
    
    parser = argparse.ArgumentParser(prog = 'Challenge_01',
                                     description = 'Program to perform the first task in the Amphora Code Challenge')
    parser.add_argument('-s', '--samples_path')
    parser.add_argument('-o', '--out_file')
    args = parser.parse_args()
    
    if args.samples_path == None:
        logger.error(f'ERROR: samples path is required!!!')
    elif args.out_file == None:
        logger.error(f'ERROR: out file is required!!!')
    else:  
        logger.info(f'listing files in folder {args.samples_path} ...')
        file_list = os.listdir(args.samples_path)
        logger.info(f'proccessing {len(file_list)} files')

        geno_table = pd.DataFrame(columns=['CHROM', 'POS', 'REF', 'ALT',])
        j = 0
        for i in file_list:
            j += 1
            if i.lower().endswith(('csv', 'vcf', 'vcf.gz')):
                if i.lower().endswith(('vcf', 'vcf.gz')):
                    format = 'vcf'
                else: 
                    format = 'csv'
                logger.info(f'({j}/{len(file_list)}) reading {format} file: {i}')
                sample_geno = import_geno(args.samples_path + i, format=format)
                sample_geno.iloc[:, -1] = geno2num(sample_geno.iloc[:, -1])
                if geno_table.empty:
                    geno_table = sample_geno
                else:
                    geno_table = geno_table.merge(sample_geno, how='left')
            else:
                logger.warning(f'({j}/{len(file_list)}) file {i} have an unknown extension')
        logger.info(f'Saving genotype table into file: {args.out_file}')
        geno_table.to_csv((args.out_file), index=False)
