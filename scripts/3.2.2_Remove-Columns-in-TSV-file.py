#!/usr/bin/python


"""
#Script to Pars Remove Columns from TSC/CSV files...

#####------Inputs-------
# ./Sript.py --Workdir
# python 4.2.2_Remove-Columns-in-TSV-file.py input.csv output.csv /home/gala0002/proj/proj_dir/NG-14833_2.1_STAR_Genome-CDS-ass/
# python 4.2.2_Remove-Columns-in-TSV-file.py Express_CDS_Pas-Reg-Tit.gene.counts.matrix_head Express_CDS_Pas-Reg-Tit.gene.counts.matrix_pars /home/gala0002/proj/proj_dir/NG-14833_2.1_STAR_Genome-CDS-ass/
#sed -i 's/_rep6/_rep4/g' Express_CDS_Pas-Reg-Tit.gene.counts.matrix_4rep

"""
import csv
import sys

input_file = sys.argv[1] #'input.csv'
output_file = sys.argv[2] #'output.csv'
WORKDIR = sys.argv[3]

cols_to_remove = [4,5,10,11,16,17] # Column indexes to be removed (starts at 0)
cols_to_remove = sorted(cols_to_remove, reverse=True) # Reverse so we remove from the end first
row_count = 0 # Current amount of rows processed

with open(WORKDIR+input_file, "r") as source:
    reader = csv.reader(source, delimiter='\t')
    with open(WORKDIR+output_file, "w") as result:
        writer = csv.writer(result, delimiter='\t')
        for row in reader:
            row_count += 1
            #print row_count
            print('\r{0}'.format(row_count)) # Print rows processed
            for col_index in cols_to_remove:
                del row[col_index]
            writer.writerow(row)

