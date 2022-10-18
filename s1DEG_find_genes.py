#-----------------Description--------------------

#-----------------Libraries----------------------
import pandas as pd

#----------------Functions-----------------------
#---------------Load Data------------------------
file = open(f'file_name.txt', "r+")
file_name = file.readline()
file.close()

DEG_shLUC = pd.read_csv(f'Arkusze/{file_name} (Ekspr+MM)/Ekspresja/DEG_{file_name}_vs_shLUC_MS.csv', delimiter= ';')
DEG_WT = pd.read_csv(f'Arkusze/{file_name} (Ekspr+MM)/Ekspresja/DEG_{file_name}_vs_WT_MS.csv', delimiter= ';')

#---------------CODE-----------------------------
#----------FIND SHARING GENES--------------------

print('-------------------------')
print('START: FIND SHARING GENES')
shLUC_genes = DEG_shLUC['Gene_Symbol'].tolist()  # get genes names
WT_genes = DEG_WT['Gene_Symbol'].tolist()

shared_genes = list()  # declare lists for genes
shared_genes = list(set(shLUC_genes) & set(WT_genes))  # get matching values between lists
shared_genes = sorted(shared_genes)
print(f'END: FIND SHARING GENES; Numb: {len(shared_genes)}')
print('-------------------------')

# -----------Find up-down regulated genes--------
print('START: FIND UP-DOWN EKSPR GENES')
DEG_ekspr = pd.DataFrame()
DEG_ekspr['Gene_Symbol'] = '-'
DEG_ekspr['Ekspresion_shLUC'] = '-'
DEG_ekspr['Ekspresion_WT'] = '-'
DEG_ekspr['Ekspresja up/down (shLuc | WT)'] = '-'

row_cnt = 0  # Copy values from DEW_shLUC/WT df to another
for x in shared_genes:
    shLUC_indx = DEG_shLUC.index[(DEG_shLUC['Gene_Symbol'] == x)].tolist()  # extract gene index from DF
    shLUC_indx = shLUC_indx[0]  # extract numerical index val
    WT_indx = DEG_WT.index[(DEG_WT['Gene_Symbol'] == x)].tolist()
    WT_indx = WT_indx[0]

    DEG_ekspr.at[row_cnt, 'Gene_Symbol'] = x
    DEG_ekspr.at[row_cnt, 'Ekspresion_shLUC'] = DEG_shLUC.at[
        shLUC_indx, 'Ekspresion_shLUC']  # copy value from one DF to another
    DEG_ekspr.at[row_cnt, 'Ekspresion_WT'] = DEG_WT.at[WT_indx, 'Ekspresion_WT']  # !!! be careful with index !!!
    row_cnt += 1  # increment row counter

# ------------find up|down reg genes------------------
DEG_eksp_row_nmb = len(DEG_ekspr.index)  # get nmb of rows
DEG_eksp_row_nmb_range = range(0, DEG_eksp_row_nmb)  # create range lst
eksp_shLUC = DEG_ekspr['Ekspresion_shLUC'].apply(lambda x: x.replace(',', '.'))  # replace ',' in every x to '.'
eksp_WT = DEG_ekspr['Ekspresion_WT'].apply(
    lambda x: x.replace(',', '.'))  # !!! remember to check how numbers are writen,

for idx, ex_Luc, exp_WT in zip(DEG_eksp_row_nmb_range, eksp_shLUC, eksp_WT):  # iterare in 3 lists
    ud_Luc = 'up' if float(ex_Luc) > 0 else 'down'  # if in one line
    ud_WT = 'up' if float(exp_WT) > 0 else 'down'
    DEG_ekspr.at[idx, 'Ekspresja up/down (shLuc | WT)'] = f'{ud_Luc} | {ud_WT}'  # insert new value to existing cell

# ----------count nmb of up/down/up-down reg genes
chck_ud = DEG_ekspr['Ekspresja up/down (shLuc | WT)'].apply(lambda x: x.replace(' ', '').split(
    '|'))  # replace+split for every val in this column, then get list of items from col
cnt_up = 0
cnt_down = 0
cnt_up_down = 0
for x in chck_ud:
    if x[0] == 'up' and x[1] == 'up':
        cnt_up += 1
    elif x[0] == 'down' and x[1] == 'down':
        cnt_down += 1
    elif (x[0] == 'up' and x[1] == 'down') or (x[0] == 'down' and x[1] == 'up'):
        cnt_up_down += 1

print(f'END: FIND UP-DOWN REG GENES; UP:{cnt_up}; DOWN:{cnt_down}; UP/DOWN:{cnt_up_down};')
print('-------------------------')
# ---------------SAVE FILE------------------------
DEG_ekspr.to_csv(f'Arkusze/{file_name} (Ekspr+MM)/Ekspresja/Results/DEG_{file_name}.csv', sep=';', header=True)
