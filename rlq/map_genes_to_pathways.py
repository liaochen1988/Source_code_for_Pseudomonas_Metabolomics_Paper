import pandas as pd
from bioservices.kegg import KEGG
import re
import numpy as np

# Read RNA-seq table
df = pd.read_csv('processed_count_table_0915.csv', index_col=0)

# Map genes to pathways
k = KEGG()
KEGG_PA_strains = ["pae", "paev", "paei", "pau",  "pap",
                   "pag", "paf",  "pnc",  "paeb", "pdk",
                   "psg", "prp",  "paep", "paer", "paem",
                   "pael","paes", "paeu", "paeg", "paec",
                   "paeo"]

dict_PAgene_KEGGpathway = {} # mappnig of PA gene name to pathway
KEGGpathway_Legend = {} # legend of KEGG pathways

count = 0
lines2append = []
PA_genes_all = dict(zip(df[['Gene','LocusTag']].drop_duplicates().Gene,df[['Gene','LocusTag']].drop_duplicates().LocusTag))
PA_genes_all = {k:v for k,v in PA_genes_all.items() if str(v)!='nan'}
for gene_id, PA14_locus in PA_genes_all.items():
    if (count % 100 == 0):
        print(count)

    gene_id_suffix_removed = re.sub(r'_\d+$', '', gene_id)
    curr_line = [gene_id, gene_id_suffix_removed, PA14_locus]
    for strain_to_map in KEGG_PA_strains:
        try:
            res = k.get_pathway_by_gene(gene_id_suffix_removed, strain_to_map)
            if (res is None) and (strain_to_map == "pau"):
                res = k.get_pathway_by_gene(PA14_locus, strain_to_map)
            if res is not None:
                curr_line.append((',').join(list(res.keys())))
                KEGGpathway_Legend = {**KEGGpathway_Legend, **res}
            else:
                curr_line.append('')
        except:
            curr_line.append('')
    lines2append.append(curr_line)
    print(curr_line)
    count += 1

df_kegg_pathway = pd.DataFrame(lines2append, index = PA_genes_all, columns = ["Gene", "Gene_wo_suffix", "PA14_locus"]+KEGG_PA_strains)
df_kegg_pathway.to_csv('kegg_pathway_mapping_raw.csv')

df_KEGGpathway_legend = pd.DataFrame.from_dict(KEGGpathway_Legend, orient='index')
df_KEGGpathway_legend.columns = ['Name']
df_KEGGpathway_legend.to_csv('KEGGpathway_Legend.csv')
