import BBH_python3
import pandas as pd
import os


subject = 'Pseudomonas_putida_KT2440_110'
query = 'Pseudomonas_aeruginosa_UCBPP_PA14'

BBH_python3.make_blast_db(subject)
BBH_python3.make_blast_db(query)
BBH_python3.get_bbh(query, subject, in_folder='bbh', outfile='')
##BBH.make_blast_db(subject)

# BBH_python3.get_bbh(query, subject, in_folder='bbh', outfile='')

##print(BBH.get_gene_lens(query, in_folder='prots'))