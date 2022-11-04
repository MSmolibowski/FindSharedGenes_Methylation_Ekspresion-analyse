#-----------------Description--------------------

#-----------------Libraries----------------------
import pandas as pd
import openpyxl
#----------------Functions-----------------------
#********************************************************************
#********************************************************************
#---------------Load Data------------------------
file = open(f'file_name.txt', "r+")         #read file name
file_name = file.readline()
file.close()

Finall_df = pd.read_csv(f'Arkusze/{file_name} (Ekspr+MM)/Results/{file_name}_final.csv', delimiter= ';')
#---------------CODE-----------------------------

#---------------SAVE FILE------------------------
Finall_df.to_excel(f'Arkusze/{file_name} (Ekspr+MM)/Results/ex{file_name}_final.xlsx', sheet_name= 'Final_df')              #save as excel file