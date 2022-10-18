#-----------------Description--------------------

#-----------------Libraries----------------------
import pandas as pd
import functools
import operator
import numpy as np
import itertools
import json

#----------------Functions-----------------------
#********************************************************************
def edit_str_from_list(edit_df, col_name):
    edit_df[col_name] = edit_df[col_name].apply(lambda x: x.replace('[', '').replace(']', '').replace("'", '').replace("''", '')) #!!! never replace ' ' (white space, sometimes list are separated by it (probably after uniqe()

    return edit_df
#********************************************************************
#---------------Load Data------------------------
file = open(f'file_name.txt', "r+")
file_name = file.readline()
file.close()

DEG_data = pd.read_csv(f'Arkusze/{file_name} (Ekspr+MM)/Ekspresja/Results/DEG_{file_name}.csv', delimiter= ';')
MM_data = pd.read_csv(f'Arkusze/{file_name} (Ekspr+MM)/Macierze Metylacji/Results/MM_{file_name}_shLUC_WT.csv', delimiter= ';')
Merge_data = pd.read_csv(f'Arkusze/{file_name} (Ekspr+MM)/Results/Merg_{file_name}.csv', delimiter= ';')
#---------------CODE-----------------------------

#1. Get info from: DEG_df
print('---------------')
print(f'START: GET INFO FROM: DEG_df')
DEG_genes_nmb = len(DEG_data.index)     #numb of genes in DEG_df
DEG_data['halp Exp'] = DEG_data['Ekspresja up/down (shLuc | WT)'].apply(lambda x: x.replace(' ','').split('|'))
DEG_data['halp Exp'] = DEG_data['halp Exp'].apply(lambda x: sorted(x))  #sort values in list

for row in range(0, DEG_genes_nmb):             #create string up/down/up_down
    exp_chck = DEG_data.at[row, 'halp Exp']
    hlp_exp = ''
    if exp_chck[0] == exp_chck[1]:
        hlp_exp = exp_chck[0]
    else:
        hlp_exp = "_".join(exp_chck)
    DEG_data.at[row, 'halp Exp'] = hlp_exp

#get uniq val from 'help Exp' and then count them
up_down_lst = DEG_data['halp Exp'].tolist()         #list of exp names
up_down_uniq = np.unique(up_down_lst)               #uniq names

DEG_answer = {'DEG_genes' : DEG_genes_nmb}          #dictionary for answer
for count_exp in up_down_uniq:                      #count uniq names
    DEG_answer[count_exp] = up_down_lst.count(count_exp)

print(f'\tDEG info: {DEG_answer}')
print(f'END: GET INFO FROM: DEG_df')
print('---------------')

#2. Get info from: MM_df
print(f'START: GET INFO FROM: MM_df')

#create met cmt list
MM_genes_nmb = len(MM_data.index)
MM_data = edit_str_from_list(MM_data, 'met_cmt_uniq')
MM_data['met_cmt_uniq'] = MM_data['met_cmt_uniq'].apply(lambda x: x.split(' '))
MM_data['met_cmt_uniq'] = MM_data['met_cmt_uniq'].apply(lambda x: sorted(x))
MM_data['hlp_met_cmt'] = MM_data['met_cmt_uniq'].apply(lambda x: '_'.join(x))
mm_met_cmd_lst = MM_data['hlp_met_cmt'].tolist()
for x in range(0, len(mm_met_cmd_lst)):        #find in list element's with first letter '_', then copy this string from his second letter
                                            #and replace in list
    chck = mm_met_cmd_lst[x]
    if chck[0] == '_':
        chck = chck[1:]
        mm_met_cmd_lst[x] = chck

met_cmd_uniq = np.unique(mm_met_cmd_lst).tolist()      #!!! remember np.uniqe().tolist() !!!, to avoid list with whate space as separator

MM_answer = {'MM_genes' : MM_genes_nmb}          #dictionary for answer
for count_met in met_cmd_uniq:                      #count uniq names
    MM_answer[count_met] = mm_met_cmd_lst.count(count_met)
print(f'\tMM info: {MM_answer}')
print(f'END: GET INFO FROM: MM_df')
print('---------------')
#3. Get info from: MM_df
print(f'START: GET INFO FROM: Merge_df')
Merge_gene_nmb = len(Merge_data.index)
Merge_data = edit_str_from_list(Merge_data, 'met_cmt_uniq')                             #Create uniq merged names of met cmmnt's
Merge_data['met_cmt_uniq'] = Merge_data['met_cmt_uniq'].apply(lambda x: x.split(' '))
Merge_data['met_cmt_uniq'] = Merge_data['met_cmt_uniq'].apply(lambda x: sorted(x))
Merge_data['met_cmt_uniq'] = Merge_data['met_cmt_uniq'].apply(lambda x: '_'.join(x))
mrg_met_cmd_lst = Merge_data['met_cmt_uniq'].tolist()

for x in range(0, len(mrg_met_cmd_lst)):
    chck = mrg_met_cmd_lst[x]
    if chck[0] == '_':
        chck = chck[1:]
        mrg_met_cmd_lst[x] = chck

Merge_data['met_cmt_uniq'] = mrg_met_cmd_lst        #insert uniq merged met cmt to column
Merge_data['Met_Exp'] = Merge_data['Ekspresja up/down (shLuc | WT)'].apply(lambda x: x.replace(' ','').split('|'))
Merge_data['Met_Exp'] = Merge_data['Met_Exp'].apply(lambda x: sorted(x))  #sort values in list

for row in range(0, Merge_gene_nmb):             #create string up/down/up_down
    chck_exp = Merge_data.at[row, 'Met_Exp']
    hlp_exp = ''
    if chck_exp[0] == chck_exp[1]:
        hlp_exp = chck_exp[0]
    else:
        hlp_exp = "_".join(chck_exp)
    Merge_data.at[row, 'Met_Exp'] = hlp_exp

Merge_data['Met_Exp'] = Merge_data['met_cmt_uniq'] + '|' + Merge_data['Met_Exp']
mrg_met_exp_lst = Merge_data['Met_Exp'].tolist()
mrg_met_exp_uniq = np.unique(mrg_met_exp_lst).tolist()
#print(mrg_met_exp_uniq)

Merge_answer = {'Merge_genes' : Merge_gene_nmb}
for count_met_exp in mrg_met_exp_uniq:
    Merge_answer[count_met_exp] = mrg_met_exp_lst.count(count_met_exp)

print(f'\tMerge info: {Merge_answer}')
print(f'END: GET INFO FROM: Merge_df')
print('---------------')

print(f'START: SAVE INFO TO TXT FILE')

with open(f'Arkusze/{file_name} (Ekspr+MM)/Results/{file_name}_info.txt', 'w') as convert_file:
    convert_file.writelines(json.dumps(DEG_answer))
    convert_file.writelines('\n')
    convert_file.writelines(json.dumps(MM_answer))
    convert_file.writelines('\n')
    convert_file.writelines(json.dumps(Merge_answer))

print(f'END: SAVE INFO TO TXT FILE')
print('---------------')

#---------------SAVE FILE------------------------
Merge_data.to_csv(f'Arkusze/{file_name} (Ekspr+MM)/Results/Merg_{file_name}.csv', sep= ';', header = True)