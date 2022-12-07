import os
import sys
import pandas as pd
from functions import *
# pd.options.mode.chained_assignment = None  # default='warn'

if __name__ == "__main__":

    samples_path = sys.argv[1]
    out_dir = sys.argv[2]

    file_list = os.listdir(samples_path)
    vcf_files = [s for s in file_list if "vcf" in s.lower()]
    csv_files = [s for s in file_list if "csv" in s.lower()]

    geno_table = pd.DataFrame(columns=["CHROM", "POS", "REF", "ALT",])
    j = 0
    for i in vcf_files:
        sample_file = samples_path + i
        j += 1
        print(j)
        sample_geno = import_geno(sample_file, format="vcf")
        # print(sample_geno)
        sample_geno.iloc[:, -1] = geno2num(sample_geno.iloc[:, -1])
        if geno_table.empty:
            geno_table = sample_geno
        else:
            geno_table = geno_table.merge(sample_geno, how="left")
    for i in csv_files:
        sample_file = samples_path + i
        j += 1
        print(j)
        sample_geno = import_geno(sample_file, format="csv")
        sample_geno.iloc[:, -1] = geno2num(sample_geno.iloc[:, -1])
        if geno_table.empty:
            geno_table = sample_geno
        else:
            geno_table = geno_table.merge(sample_geno, how="left")
    geno_table.to_csv((out_dir + "geno_table.csv"), index=False)
    print(geno_table.shape)


