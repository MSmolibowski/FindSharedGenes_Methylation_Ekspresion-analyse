#-----------------Description--------------------

#-----------------Libraries----------------------
import pandas as pd
import functools
import operator
import numpy as np

#----------------Functions-----------------------
#********************************************************************
def create_gene_lst(edit_df):
    gene_lst = functools.reduce(operator.iadd, edit_df['Gene_Symbol_uniq_help'].tolist(), [])  # get list of list from ['Localization'] col and merge all to one list
    # also we can use operator.iconcat than operator.iadd, both equal fast with many small list, or few long lists
    gene_lst = np.unique(gene_lst)
    return gene_lst

def del_unnamed(edit_df):
    del edit_df['Unnamed: 0']
    return edit_df

def create_string(edit_df):
    edit_df['Gene_Symbol_uniq'] = edit_df['Gene_Symbol_uniq'].apply(lambda x: x.replace('[','').replace(']','').replace("'",''))
    #print(type(edit_df.at[0, 'Gene_Symbol_uniq']))
    return edit_df

def reindex_df(edit_df):
    edit_df.to_csv(f'Arkusze/{file_name} (Ekspr+MM)/Macierze Metylacji/reindex/DF_reindex.csv', sep=';', header=True)
    edit_df = pd.read_csv(f'Arkusze/{file_name} (Ekspr+MM)/Macierze Metylacji/reindex/DF_reindex.csv', delimiter=';')
    del edit_df['Unnamed: 0']
    return edit_df

def find_and_fill(new_df, shLUC_df, WT_df):
    row_nmb = len(new_df.index)

    for row in range(0, row_nmb):
        gene = new_df.at[row, 'Gene_Symbol']                                            #get gene shared gene name (gene in shLUC and WT)
        shLUC_indx = shLUC_df.index[(shLUC_df['Gene_Symbol_uniq'] == gene)].tolist()    #get index of row with matchng gene name in shLUC_df
        WT_indx = WT_df.index[(WT_df['Gene_Symbol_uniq'] == gene)].tolist()             #get index of row with matchng gene name in WT_df

        LUC_probe = list()  #list for values
        LUC_loc = list()
        LUC_met = list()
        WT_probe = list()
        WT_loc = list()
        WT_met = list()
        met_cmmt =list()

        if len(shLUC_indx) > 0:
            for indx_L in shLUC_indx:                                   #for every index, find and add values to list
                LUC_probe.append(shLUC_df.at[indx_L, 'TargetID'])
                LUC_loc.append(shLUC_df.at[indx_L, 'Localization'])
                LUC_met.append(shLUC_df.at[indx_L, 'Methylation_val'])
                met_cmmt.append(shLUC_df.at[indx_L, 'Methylation_com'])

            new_df.at[row, 'probesID_shLUC'] = LUC_probe                #add list to cells
            new_df.at[row, 'Localization_shLUC'] = LUC_loc
            new_df.at[row, 'met_val_shLUC'] = LUC_met

        if len(WT_indx) > 0:
            for indx_W in WT_indx:
                WT_probe.append(WT_df.at[indx_W, 'TargetID'])
                WT_loc.append(WT_df.at[indx_W, 'Localization'])
                WT_met.append(WT_df.at[indx_W, 'Methylation_val'])
                met_cmmt.append(WT_df.at[indx_W, 'Methylation_com'])

            new_df.at[row, 'probesID_WT'] = WT_probe
            new_df.at[row, 'Localization_WT'] = WT_loc
            new_df.at[row, 'met_val_WT'] = WT_met
            new_df.at[row, 'met_cmt'] = met_cmmt
        #print(f'{gene}; shLUC: {shLUC_indx}; WT {WT_indx}')
    return new_df

def get_rest_val(final_df, edit_df, switcher):
    prob = 'probesID_' + switcher
    local = 'Localization_' + switcher
    met = 'met_val_' + switcher
    # =====loop to get gene from row and find proper index in Finall DF
    help_len = len(edit_df.index)
    edit_df['Gene_Symbol_uniq'] = edit_df['Gene_Symbol_uniq'].apply(lambda x: x.replace('[','').replace(']', '').replace("'", '').split(' '))

    for row in range(0, help_len):  #
        gene = edit_df.at[row, 'Gene_Symbol_uniq']  # list of genes
        for g in gene:
            merg_indx = final_df.index[(final_df['Gene_Symbol'] == g)].tolist()
            probeID_lst = list()
            loc_lst = list()
            mev_val_lst = list()
            met_cmmt = list()

            if len(merg_indx) > 0:
                probeID_lst.append(final_df.at[merg_indx[0], prob]) #probeID list from Merge_DF
                probeID_add = edit_df.at[row, 'TargetID']           #value to add
                probeID_lst.append(probeID_add)                     #add value to list
                final_df.at[merg_indx[0], prob]  = probeID_lst   #insert new val

                # #Localization
                loc_lst.append(final_df.at[merg_indx[0], local])  # probeID list from Merge
                loc_add = edit_df.at[row, 'Localization']
                loc_lst.append(loc_add)
                final_df.at[merg_indx[0], local] = loc_lst  # insert new val:

                #met_value
                mev_val_lst.append(final_df.at[merg_indx[0], met])  # probeID list from Merge
                met_val_add = edit_df.at[row, 'Methylation_val']
                mev_val_lst.append(met_val_add)
                final_df.at[merg_indx[0], met] = mev_val_lst  # insert new val:

                #met_cmmnt
                met_cmmt.append(final_df.at[merg_indx[0], 'met_cmt'])  # probeID list from Merge
                met_cmmt_add = edit_df.at[row, 'Methylation_com']
                met_cmmt.append(met_cmmt_add)
                final_df.at[merg_indx[0], 'met_cmt'] = met_cmmt  # insert new val:

    return final_df
#********************************************************************
#---------------Load Data------------------------
file = open(f'file_name.txt', "r+")
file_name = file.readline()
file.close()

MM_shLUC = pd.read_csv(f'Arkusze/{file_name} (Ekspr+MM)/Macierze Metylacji/Results/MM_shLUC_{file_name}_shared_only.csv', delimiter= ';')
MM_WT = pd.read_csv(f'Arkusze/{file_name} (Ekspr+MM)/Macierze Metylacji/Results/MM_WT_{file_name}_shared_only.csv', delimiter= ';')
MM_shLUC = del_unnamed(MM_shLUC)
MM_WT = del_unnamed(MM_WT)
#---------------CODE-----------------------------
#1. Edit 'Gene_Symbol_uniq' column to list
#create help col
MM_shLUC['Gene_Symbol_uniq_help'] = MM_shLUC['Gene_Symbol_uniq'].apply(lambda x: x.replace('[','').replace(']', '').replace("'", '').split(' ')) #create help col with gene list
MM_WT['Gene_Symbol_uniq_help'] = MM_WT['Gene_Symbol_uniq'].apply(lambda x: x.replace('[','').replace(']', '').replace("'", '').split(' '))

#2. Find genes shared between MM_DF
print(f'--------------------')
print(f'START: FIND SHARED GENES BETWEEN MM_DFs')
shLUC_genes = create_gene_lst(MM_shLUC)
WT_genes = create_gene_lst(MM_WT)
shared_genes = list(set(shLUC_genes)&set(WT_genes))
shared_genes = sorted(shared_genes)
print(f'MM shared genes: {len(shared_genes)}; MM_shLUC: {len(shLUC_genes)}; MM_WT: {len(WT_genes)}')
print(f'END: FIND SHARED GENES BETWEEN MM_DFs')
print(f'--------------------')

#3. create DF with uniq_gene == 1 and >1
shLUC_1 = MM_shLUC[MM_shLUC['Genes_uniq_len']==1]       #helpd df
shLUC_more = MM_shLUC[MM_shLUC['Genes_uniq_len']>1]
shLUC_1 = reindex_df(shLUC_1)                           #reindex help df
shLUC_more = reindex_df(shLUC_more)

WT_1 = MM_WT[MM_WT['Genes_uniq_len']==1]
WT_more = MM_WT[MM_WT['Genes_uniq_len']>1]
WT_1 = reindex_df(WT_1)
WT_more = reindex_df(WT_more)
#create string in 'Gene_Symbol_uniq' column
shLUC_1 = create_string(shLUC_1)
WT_1 = create_string(WT_1)

#4. Create new DF for shared genes between shLUC_MM
print(f'START: CREATE NEW DF AND FILL WITH VALUES')
MM_Merge = pd.DataFrame()
MM_Merge['Gene_Symbol'] = shared_genes
MM_Merge['probesID_shLUC'] = ''
MM_Merge['probesID_WT'] = ''
MM_Merge['Localization_shLUC'] = ''
MM_Merge['Localization_WT'] = ''
MM_Merge['met_val_shLUC'] = ''
MM_Merge['met_val_WT'] = ''
MM_Merge['met_cmt'] = ''

print(f'\t START: USE find_and_fill() function')
MM_Merge = find_and_fill(MM_Merge, shLUC_1, WT_1)
print(f'\t END: USE find_and_fill() function')

print(f'\t\t START: USE get_rest_val() function for  shLUC_df')
MM_Merge = get_rest_val(MM_Merge, shLUC_more, 'shLUC')
print(f'\t END: USE get_rest_val() function')

print(f'\t\t START: USE get_rest_val() function for  WT_df')
MM_Merge = get_rest_val(MM_Merge, WT_more, 'WT')
print(f'\t END: USE get_rest_val() function')

print(f'END: CREATE NEW DF AND FILL WITH VALUES')
print(f'--------------------')

#4. Create col with uniq names of hypo, hyper names
#print(MM_Merge['met_cmt'])
MM_Merge = reindex_df(MM_Merge)
MM_Merge['met_cmt'] = MM_Merge['met_cmt'].apply(lambda x: x.replace('[','').replace(']','').replace("'",'').replace("''",'').replace(' ',''))
MM_Merge['met_cmt'] = MM_Merge['met_cmt'].apply(lambda x: x.split(','))
MM_Merge['met_cmt_uniq'] = MM_Merge['met_cmt'].apply(lambda x: np.unique(x))


#---------------SAVE FILE------------------------
MM_Merge.to_csv(f'Arkusze/{file_name} (Ekspr+MM)/Macierze Metylacji/Results/MM_{file_name}_shLUC_WT.csv', sep = ';', header=True)