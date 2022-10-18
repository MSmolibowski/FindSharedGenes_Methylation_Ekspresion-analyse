#-----------------Description--------------------

#-----------------Libraries----------------------
import pandas as pd
import functools
import numpy as np
import operator
import numpy as np

#----------------Functions-----------------------
#********************************************************************
def find_row_without_TSS(edit_df):
    edit_df_row_nmb = len(edit_df.index)
    edit_df_TSS_indx = list()  # get rows indx with TSS1500 or TSS200
    edit_df_del = list()  # get row index without TSS1500 or TSS200

    for row in range(0, edit_df_row_nmb):
        chck_loc = edit_df.at[row, 'Localization']
        if 'TSS1500' in chck_loc:
            edit_df_TSS_indx.append(row)
        else:
            if 'TSS200' in chck_loc:
                edit_df_TSS_indx.append(row)
            else:
                edit_df_del.append(row)
    return edit_df_del
#********************************************************************
def del_non_TSS(edit_df):
    edit_df['Gene_Symbol'] = edit_df['Gene_Symbol'].apply(lambda x: x.split(';'))  # create list from string by split string by ';'
    edit_df['Localization'] = edit_df['Localization'].apply(lambda x: x.split(';'))

    loc_merg_lst = functools.reduce(operator.iadd, edit_df['Localization'].tolist(), [])  # get list of list from ['Localization'] col and merge all to one list
    # also we can use operator.iconcat than operator.iadd, both equal fast with many small list, or few long lists
    loc_merg_lst = np.unique(loc_merg_lst)  # get uniq values from list, now it is numpy array !!!
    loc_del = [x for x in loc_merg_lst if (x != 'TSS1500' and x != 'TSS200')]  # get values to del from loc (values different than TSS1500 and TSS200)

    row_nmb = len(edit_df.index)
    for row in range(0, row_nmb):
        loc_val = edit_df.at[row, 'Localization']
        gne_nme = edit_df.at[row, 'Gene_Symbol']
        check = any(item in loc_del for item in loc_val)

        if check == True:
            for del_val in loc_del:  # take val_name from list of values that have to be removed
                if del_val in loc_val:  # check if this val_name is in list from cell MM_shLuc.at[row, 'Localization']
                    del_indx = [indx for indx, val in enumerate(loc_val) if val in del_val]  # get index of matching value in list

                    # print(f'{row} ----- {del_val} ----- {loc_val}-----{del_indx}')
                    # print(f'Before del: {loc_val}')
                    loc_val = np.delete(loc_val, del_indx)  # delete items by index from localization list
                    gne_nme = np.delete(gne_nme, del_indx)  # delete items by index from gne_name list
                    # print(f'Before del: {loc_val}')
                    # print('-----------------------')
            # --------insetr new loc - gen to proper cells
            edit_df.at[row, 'Localization'] = loc_val
            edit_df.at[row, 'Gene_Symbol'] = gne_nme
    return edit_df
#********************************************************************
def del_non_sharing_genes(edit_df, sharing_gne_lst):
    row_nmb = len(edit_df.index)
    for row in range(0, row_nmb):
        edit_gen = edit_df.at[row, 'Gene_Symbol']
        edit_loc = edit_df.at[row, 'Localization']
        if len(edit_gen) > 1:
            for del_gen in edit_gen:
                if del_gen not in sharing_gne_lst:
                    gen_indx = [indx for indx, val in enumerate(edit_gen) if val in del_gen]  # get index of matching value in list
                                                    #enumerate created indx:val' if val in del_gen give me indx
                    edit_loc = np.delete(edit_loc, gen_indx)  # delete items by index from localization list
                    edit_gen = np.delete(edit_gen, gen_indx)  # delete items by index from gne_name list
            edit_df.at[row, 'Localization'] = edit_loc
            edit_df.at[row, 'Gene_Symbol'] = edit_gen

    edit_df['Gene_Symbol_uniq'] = edit_df['Gene_Symbol'].apply(lambda x: np.unique(x))
    edit_df['Genes_uniq_len'] = edit_df['Gene_Symbol_uniq'].apply(lambda x: len(x))
    edit_df.drop(edit_df[(edit_df['Genes_uniq_len'] == 0)].index, inplace=True)
    return edit_df
#********************************************************************
def genes_uniq(edit_df):
    edit_genes = functools.reduce(operator.iconcat, edit_df['Gene_Symbol'].tolist(), []) # get list of list from ['Gene_Symbol'] col and merge all to one list
    edit_df_uniq_gne = np.unique(edit_genes)  # get uniq values from list, now it is numpy array !!!
    return edit_df_uniq_gne

def reindex_df(edit_df):
    edit_df.to_csv(f'Arkusze/{file_name} (Ekspr+MM)/Macierze Metylacji/reindex/DF_reindex.csv', sep=';', header=True)
    edit_df = pd.read_csv(f'Arkusze/{file_name} (Ekspr+MM)/Macierze Metylacji/reindex/DF_reindex.csv', delimiter=';')
    del edit_df['Unnamed: 0']
    return edit_df

#---------------Load Data------------------------
file = open(f'file_name.txt', "r+")
file_name = file.readline()
file.close()

MM_shLuc = pd.read_csv(f'Arkusze/{file_name} (Ekspr+MM)/Macierze Metylacji/MM_{file_name}_vs_shLUC_MS.csv', delimiter= ';')
MM_shLuc = MM_shLuc[~MM_shLuc['Gene_Symbol'].isnull()]        #copy all rows with values in col Genes
MM_shLuc = MM_shLuc[~MM_shLuc['Localization'].isnull()] #copy all rows with values in col Genes

MM_WT = pd.read_csv(f'Arkusze/{file_name} (Ekspr+MM)/Macierze Metylacji/MM_{file_name}_vs_WT_MS.csv', delimiter= ';')
MM_WT = MM_WT[~MM_WT['Gene_Symbol'].isnull()]
MM_WT = MM_WT[~MM_WT['Localization'].isnull()]
#---------------CODE-----------------------------
#1. Find index of no TSS row's in DF
print(f'--------------------------------------')
print(f'START: FIND NON TSS1500 AND TSS200 INDEX')

shLUC_del = find_row_without_TSS(MM_shLuc)
WT_del = find_row_without_TSS(MM_WT)

print(f'END: FIND NON TSS1500 AND TSS200 INDEX')
print(f'--------------------------------------')

#2. DEL NON TSS1500 OR TSS200 ROWS
print(f'START: DEL NON TSS ROWS')
MM_shLuc.drop(shLUC_del, axis = 0, inplace= True)
MM_WT.drop(WT_del, axis = 0, inplace= True )
print(f'END: DEL NON TSS ROWS')
print(f'--------------------------------------')

# 2.1 Reindex DF
MM_shLuc = reindex_df(MM_shLuc)
MM_WT = reindex_df(MM_WT)

#3.DELETE GENES WITHOUT TSS LOC
print('START: DEL GENES WITH LOC != TSS1500 AND TSS200')
print(f'\tSTART: USE FUNCTION: del_non_TSS for MM_shLUC')
MM_shLuc = del_non_TSS(MM_shLuc)
print(f'\tEND: USE FUNCTION: del_non_TSS for MM_shLUC')
print()
print(f'\tSTART: USE FUNCTION: del_non_TSS for MM_WT')
MM_WT = del_non_TSS(MM_WT)
print(f'\tEND: USE FUNCTION: del_non_TSS for MM_WT')
print('END: DEL GENES WITH LOC != TSS1500 AND TSS200')
print(f'--------------------------------------')

#===SAVE only TSS genes:
#MM_shLuc.to_csv('shLUC_TSS_only.csv', sep=';', header=True)
#MM_WT.to_csv('shLUC_TSS_only.csv', sep=';', header=True)

#4. Find genes shared between two MM_DF's
shLUC_uniq_gne = genes_uniq(MM_shLuc)           #use function to get genes from shLUC
WT_genes_uniq_gne = genes_uniq(MM_WT)
shared_genes = list(set(shLUC_uniq_gne) & set(WT_genes_uniq_gne)) #get genes shared between shLUC-WT
shared_genes = sorted(shared_genes) #Sort gene names

#print(f'MM shared genes: {len(shared_genes)}; MM_shLUC: {len(shLUC_uniq_gne)}; MM_WT: {len(WT_genes_uniq_gne)}')
#4. DEL NON SHARING GENES IN COLUMN 'Gene_Symbol' with their location
#===============EDIT GENE COLUMNS=====================
print('START: DEL NON SHARING GENES')
print(f'\tSTART: USE FUNCTION: del_non_sharing_genes for MM_shLUC')
MM_shLuc = del_non_sharing_genes(MM_shLuc, shared_genes)
print(f'\tEND: USE FUNCTION: del_non_sharing_genes for MM_shLUC')
print()
print(f'\tSTART: USE FUNCTION: del_non_sharing_genes for MM_WT')
MM_WT = del_non_sharing_genes(MM_WT, shared_genes)
print(f'\tEND: USE FUNCTION: del_non_sharing_genes for MM_WT')
print('END: DEL NON SHARING GENES')
print(f'--------------------------------------')

# # 4.1 Reindex DF !!!check values in columns
# MM_shLuc = reindex_df(MM_shLuc)
# MM_WT = reindex_df(MM_WT)

shLUC_uniq_gne = genes_uniq(MM_shLuc)           #use function to get genes from shLUC
WT_genes_uniq_gne = genes_uniq(MM_WT)
shared_genes = list(set(shLUC_uniq_gne) & set(WT_genes_uniq_gne)) #get genes shared between shLUC-WT
shared_genes = sorted(shared_genes) #Sort gene names

print(f'MM shared genes: {len(shared_genes)}; MM_shLUC: {len(shLUC_uniq_gne)}; MM_WT: {len(WT_genes_uniq_gne)}')


#---------------SAVE FILE------------------------
MM_shLuc.to_csv(f'Arkusze/{file_name} (Ekspr+MM)/Macierze Metylacji/Results/MM_shLUC_{file_name}_shared_only.csv', sep=';', header=True)
MM_WT.to_csv(f'Arkusze/{file_name} (Ekspr+MM)/Macierze Metylacji/Results/MM_WT_{file_name}_shared_only.csv', sep=';', header=True)