#-----------------Description--------------------

#-----------------Libraries----------------------
import pandas as pd
from openpyxl import Workbook
#----------------Functions-----------------------
#********************************************************************
def str_to_lst(edit_df, col_nmes):                                         #edit string in columns to list
    for col in col_nmes:
        edit_df[col] = edit_df[col].apply(lambda x: x.replace('[', '').replace(']', '').replace("'",'').replace(' ','').split(','))
        edit_df[col] = edit_df[col].apply(lambda x: '_'.join(x))

        col_str = edit_df[col].tolist()
        for x in range(0, len(col_str)):    #delete '' from list
            chck = col_str[x]
            if chck[0] == '_':
                chck = chck[1:]
                col_str[x] = chck
        edit_df[col] = col_str
        edit_df[col] = edit_df[col].apply(lambda x: x.split('_'))

    return edit_df

def cnt_hypo_hyper(met_val):        #count number of hypo and hyper probes
    hypo_cnt = 0
    hyper_cnt = 0

    for val in met_val:
        if float(val) < 0:
            hypo_cnt+=1
        if float(val) > 0:
            hyper_cnt+=1

    answer = [hypo_cnt, hyper_cnt]
    return answer

#********************************************************************
#---------------Load Data------------------------
file = open(f'file_name.txt', "r+")         #read file name
file_name = file.readline()
file.close()

Merge_data = pd.read_csv(f'Arkusze/{file_name} (Ekspr+MM)/Results/Merg_{file_name}.csv', delimiter = ';')
del Merge_data['Unnamed: 0.1']
del Merge_data['Unnamed: 0']
DEG_data_WT = pd.read_csv(f'Arkusze/{file_name} (Ekspr+MM)/Ekspresja/DEG_{file_name}_vs_WT_MS.csv', delimiter = ';')

#---------------CODE-----------------------------
#1. CREATE FINAL DATAFRAME
Finall_df = pd.DataFrame()                          #create dinall DF
Finall_df['Gene_Symbol'] = Merge_data['Gene_Symbol'] #copy genes names
Finall_df.insert(0, 'Methyl | Expr', '-')           #add column to specific position and fill rows with '-'
Finall_df['Gene_ID'] = '-'                          #create new column
Finall_df['Transcript_ID'] = '-'
Finall_df['Description'] = '-'
Finall_df['gene_biotype'] = '-'
Finall_df['Ekspresja up/down (shLuc | WT)'] = Merge_data['Ekspresja up/down (shLuc | WT)']
Finall_df['Ekspresion_shLUC'] = Merge_data['Ekspresion_shLUC']
Finall_df['Ekspresion_WT'] = Merge_data['Ekspresion_WT']
Finall_df['Total number of probes (shLuc | WT)'] = '-'
Finall_df['probesID_shLUC'] = Merge_data['probesID_shLUC']
Finall_df['Localization_shLUC'] = Merge_data['Localization_shLUC']
Finall_df['probesID_WT'] = Merge_data['probesID_WT']
Finall_df['Localization_WT'] = Merge_data['Localization_WT']
Finall_df['Hypo probes (shLuc | WT)'] = '-'
Finall_df['Hyper probes (shLuc | WT)'] = '-'
Finall_df['Hypo : Hyper (shLuc | WT)'] = '-'
Finall_df['met_val_shLUC'] = Merge_data['met_val_shLUC']
Finall_df['met_val_WT'] = Merge_data['met_val_WT']

edit_col = ['probesID_shLUC', 'probesID_WT', 'met_val_shLUC', 'met_val_WT']
Merge_data = str_to_lst(Merge_data, edit_col)      #use function to edit string to list in columns

row_nmb = len(Finall_df.index)
for row in range(0, row_nmb):                   #loop to copy values
    gene = Finall_df.at[row, 'Gene_Symbol']     #get gene name
    merg_gne_indx = Merge_data.index[(Merge_data['Gene_Symbol'] == gene)].tolist()  #find index of gene in Merge_df
    merg_gne_indx = merg_gne_indx[0]
    deg_gne_indx = DEG_data_WT.index[(DEG_data_WT['Gene_Symbol'] == gene)].tolist() #find index of gene in DEG_df
    deg_gne_indx = deg_gne_indx[0]

    #variables for counting probes
    LUC_prob = Merge_data.at[merg_gne_indx, 'probesID_shLUC']
    WT_probe = Merge_data.at[merg_gne_indx, 'probesID_WT']
    LUC_hypo_hyper = cnt_hypo_hyper(Merge_data.at[merg_gne_indx, 'met_val_shLUC'])
    WT_hypo_hyper = cnt_hypo_hyper(Merge_data.at[merg_gne_indx, 'met_val_WT'])

    #Copy info from Merge_df
    Finall_df.at[row, 'Methyl | Expr'] = Merge_data.at[merg_gne_indx, 'Met_Exp']
    Finall_df.at[row, 'Total number of probes (shLuc | WT)'] = f'{len(LUC_prob)} | {len(WT_probe)}'
    Finall_df.at[row, 'Hypo probes (shLuc | WT)'] = f'{LUC_hypo_hyper[0]} | {WT_hypo_hyper[0]}'
    Finall_df.at[row, 'Hyper probes (shLuc | WT)'] = f'{LUC_hypo_hyper[1]} | {WT_hypo_hyper[1]}'
    Finall_df.at[row, 'Hypo : Hyper (shLuc | WT)'] = f'shLUC: {LUC_hypo_hyper[0]}:{LUC_hypo_hyper[1]} | WT: {WT_hypo_hyper[0]}:{WT_hypo_hyper[1]}'

    #Copy info from DEG_df
    Finall_df.at[row,'Gene_ID'] = DEG_data_WT.at[deg_gne_indx,'Gene_ID']
    Finall_df.at[row,'Transcript_ID'] = DEG_data_WT.at[deg_gne_indx,'Transcript_ID']
    Finall_df.at[row,'Description'] = DEG_data_WT.at[deg_gne_indx,'Description']
    Finall_df.at[row,'gene_biotype'] = DEG_data_WT.at[deg_gne_indx,'gene_biotype']

Finall_df['Ekspresion_shLUC'] = Merge_data['Ekspresion_shLUC'].apply(lambda x: str(x))      #change value in column to string (if not values will be saved as 4,6 not 4.6 !!!
Finall_df['Ekspresion_WT'] = Merge_data['Ekspresion_WT'].apply(lambda x: str(x))

save_path_csv = f'Arkusze/{file_name} (Ekspr+MM)/Results/{file_name}_final.csv'
save_path_xlms = f'Arkusze/{file_name} (Ekspr+MM)/Results/ex{file_name}_final.xlsx'
#---------------SAVE FILE------------------------
Finall_df.to_csv(save_path_csv, sep=';', header=True) #save as .csv
Finall_df.to_excel(save_path_xlms)              #save as excel file
print(f'Created Final .csv file named: {file_name}_final and excel file named: ex{file_name}_final')