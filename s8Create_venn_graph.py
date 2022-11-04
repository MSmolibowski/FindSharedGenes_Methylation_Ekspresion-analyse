#-----------------Description--------------------

#-----------------Libraries----------------------
import pdb
from matplotlib import pyplot as plt
from venn import venn

import pandas as pd
from matplotlib_venn import venn2, venn2_circles, venn2_unweighted
from matplotlib_venn import venn3, venn3_circles

#----------------Functions-----------------------
#********************************************************************
#********************************************************************
#---------------Load Data------------------------
file_name = input("Type file name ex. 'SKMES_sh714': ")

DEG_data = pd.read_csv(f'Arkusze/{file_name} (Ekspr+MM)/Ekspresja/Results/DEG_{file_name}.csv', delimiter= ';')
del DEG_data['Unnamed: 0']
MM_data = pd.read_csv(f'Arkusze/{file_name} (Ekspr+MM)/Macierze Metylacji/Results/MM_{file_name}_shLUC_WT.csv', delimiter= ';')
del MM_data['Unnamed: 0']

#sets of genes //find witch genes are down/up/hyper/hypo and add them to set's
DOWN = set()
UP = set()
HYPER = set()
HYPO = set()

DEG_row = len(DEG_data.index)
for row in range(0, DEG_row):
    gene = DEG_data.at[row, 'Gene_Symbol']
    if "up" in DEG_data.at[row, 'Ekspresja up/down (shLuc | WT)']:
        UP.add(gene)
    if "down" in DEG_data.at[row, 'Ekspresja up/down (shLuc | WT)']:
        DOWN.add(gene)

MM_row = len(MM_data.index)
for row in range(0, MM_row):
    gene = MM_data.at[row, 'Gene_Symbol']
    if "hypo" in MM_data.at[row, 'met_cmt_uniq']:
        HYPO.add(gene)
    if "hyper" in MM_data.at[row, 'met_cmt_uniq']:
        HYPER.add(gene)

data_dict = {                   #create dictionary with sets
    "High expression" : UP,
    "Low expression" : DOWN,
    "Hypermethylation" : HYPER,
    "Hypomethylation" : HYPO,
}

#CREATE VENN GRAPH
plt.figure()                    #start creating figure
venn(data_dict)                 #CREATE VENN GRAPH
plt.title(f'{file_name}')
plt.savefig(f'Arkusze/{file_name} (Ekspr+MM)/Results/{file_name}_venn.png') #save venn graph



