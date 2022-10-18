# -----------------Description--------------------

# -----------------Libraries----------------------
import pandas as pd
import functools
import operator
import numpy as np
import itertools

#pd.options.mode.chained_assignment = None  # default='warn'
# ----------------Functions-----------------------
# ********************************************************************
def genes_uniq(edit_df):
    edit_genes = functools.reduce(operator.iconcat, edit_df['Gene_Symbol'].tolist(), [])  # get list of list from ['Gene_Symbol'] col and merge all to one list
    edit_df_uniq_gne = np.unique(edit_genes)  # get uniq values from list, now it is numpy array !!!
    return edit_df_uniq_gne

def edit_str_to_lst(edit_df):
    edit_df['Gene_Symbol'] = edit_df['Gene_Symbol'].apply(lambda x: x.replace('[', '').replace(']', '').replace("'", '').replace(' ', ''))
    edit_df['Gene_Symbol'] = edit_df['Gene_Symbol'].apply(lambda x: x.split(','))
    edit_df['Gene_Symbol_uniq'] = edit_df['Gene_Symbol_uniq'].apply(lambda x: x.replace('[', '').replace(']', '').replace("'", ''))

    return edit_df

# ********************************************************************
# ---------------Load Data------------------------
file = open(f'file_name.txt', "r+")
file_name = file.readline()
file.close()

DEG_data = pd.read_csv(f'Arkusze/{file_name} (Ekspr+MM)/Ekspresja/Results/DEG_{file_name}.csv', delimiter=';')
MM_data = pd.read_csv(f'Arkusze/{file_name} (Ekspr+MM)/Macierze Metylacji/Results/MM_{file_name}_shLUC_WT.csv', delimiter=';')

# ---------------CODE-----------------------------
#1. FIND GENES SHARED BETWEEN MM-DEG
print(f'--------------------')
print(f'START: FIND SHARED GENES AND CREATE NEW DF')

DEG_genes = DEG_data['Gene_Symbol'].tolist()  # get DEG genes
MM_gene = MM_data['Gene_Symbol'].tolist()

shared_genes = list(set(DEG_genes) & set(MM_gene))  # get shared genes between DEG and MM
shared_genes = sorted(shared_genes)  # sort genes names
#2. Create DF for DEG_MM
DEG_MM_Merge = pd.DataFrame()
DEG_MM_Merge['Gene_Symbol'] = shared_genes
DEG_MM_Merge['probesID_shLUC'] = ''
DEG_MM_Merge['probesID_WT'] = ''
DEG_MM_Merge['Localization_shLUC'] = ''
DEG_MM_Merge['Localization_WT'] = ''
DEG_MM_Merge['met_val_shLUC'] = ''
DEG_MM_Merge['met_val_WT'] = ''
DEG_MM_Merge['met_cmt'] = ''
DEG_MM_Merge['met_cmt_uniq'] = ''
DEG_MM_Merge['Ekspresion_shLUC'] = ''
DEG_MM_Merge['Ekspresion_WT'] = ''
DEG_MM_Merge['Ekspresja up/down (shLuc | WT)'] = ''

print(f'END: FIND SHARED GENES AND CREATE NEW DF')
print(f'--------------------')
print(f'START: FIND AND COPY VALUES')
row_nmb = len(DEG_MM_Merge.index)

for row in range(0, row_nmb):
    gene = DEG_MM_Merge.at[row, 'Gene_Symbol']
    MM_gene_indx = MM_data.index[(MM_data['Gene_Symbol'] == gene)].tolist()
    MM_gene_indx = MM_gene_indx[0]

    DEG_gne_indx = DEG_data.index[(DEG_data['Gene_Symbol'] == gene)].tolist()
    DEG_gne_indx = DEG_gne_indx[0]
    #print(f'Gene {gene}; DEG: {DEG_gne_indx}; MM {MM_gene_indx}')
    #MM info
    DEG_MM_Merge.at[row,'probesID_shLUC'] = MM_data.at[MM_gene_indx, 'probesID_shLUC']
    DEG_MM_Merge.at[row,'probesID_WT'] = MM_data.at[MM_gene_indx, 'probesID_WT']
    DEG_MM_Merge.at[row,'Localization_shLUC'] = MM_data.at[MM_gene_indx, 'Localization_shLUC']
    DEG_MM_Merge.at[row,'Localization_WT'] = MM_data.at[MM_gene_indx, 'Localization_WT']
    DEG_MM_Merge.at[row,'met_val_shLUC'] = MM_data.at[MM_gene_indx, 'met_val_shLUC']
    DEG_MM_Merge.at[row,'met_val_WT'] = MM_data.at[MM_gene_indx, 'met_val_WT']
    DEG_MM_Merge.at[row,'met_cmt'] = MM_data.at[MM_gene_indx, 'met_cmt']
    DEG_MM_Merge.at[row, 'met_cmt_uniq'] = MM_data.at[MM_gene_indx, 'met_cmt_uniq']

    #DEG info
    DEG_MM_Merge.at[row,'Ekspresion_shLUC'] = DEG_data.at[DEG_gne_indx,'Ekspresion_shLUC']
    DEG_MM_Merge.at[row,'Ekspresion_WT'] = DEG_data.at[DEG_gne_indx,'Ekspresion_WT']
    DEG_MM_Merge.at[row,'Ekspresja up/down (shLuc | WT)'] = DEG_data.at[DEG_gne_indx,'Ekspresja up/down (shLuc | WT)']

#do przemyslenia rozwiazanie: usunięcie wierszy z genami niewspolnymi a następnie połączenie obu arkuszy razem

print(f'END: FIND AND COPY VALUES')
print(f'--------------------')
#---------------SAVE FILE------------------------
DEG_MM_Merge.to_csv(f'Arkusze/{file_name} (Ekspr+MM)/Results/Merg_{file_name}.csv', sep = ';', header = True)
