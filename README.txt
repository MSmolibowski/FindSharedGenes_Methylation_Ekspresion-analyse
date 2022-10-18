DESCRIPTOIN:
- find which genes are shared between 2 datsheet's from RNA-Seq (DEG_df), and 2 Methylation Matrix (MM_df);
- find information how much genes are up/down regulated (DEG_df),
- find information how much genes are hyper/hypo methylated (MM_df),
- find information how much genes are hypo-up/hypo-down/hyper-up etc. (DEG_df + MM_df = Merg_df)
------------------------

PREPARE FOLDERS:
- create folders in the same place where scripts
  Arkusze/{file_name} (Ekspr+MM)/Ekspresja/Results
  Arkusze/{file_name} (Ekspr+MM)/Macierze Metylacji/reindex
  Arkusze/{file_name} (Ekspr+MM)/Macierze Metylacji/Results
  Arkusze/{file_name} (Ekspr+MM)/Results
------------------------

PREPARE DATA:
- create .csv files separated by ';',
- name your DEG_df files: DEG_{file_name}_vs_shLUC_MS, DEG_{file_name}_vs_WT_MS,
- name your MM_df files: MM_{file_name}_vs_shLUC_MS, MM_{file_name}_vs_WT_MS, (columns: TargetID, Methylation_val, Gene_Symbol, Localization)
- place DEG_df .csv file separated by ';' in in 'Arkusze/{file_name} (Ekspr+MM)/Ekspresja'
- place MM_df .csv file separated by ';' in in 'Arkusze/{file_name} (Ekspr+MM)/Macierze Metylacji'

SCRIPTS DESCRIPTION:
*s0Main.py*
~ take file name from user example: SKMES_sh714
~ save it to .txt file,
~ run scripts,
------------------------

*s1DEG_find_genes*
~ find shared genes between DEG files,
~ create new .csv file named 'DEG_{file_name}' with shared genes in 'Arkusze/{file_name} (Ekspr+MM)/Macierze Metylacji/Results',
------------------------

*s2MM_filter*
~ find genes NOT shared between MM_ files,
~ delete rows with not shared genes from MM,
~ find genes without localization in TSS1500 and TSS200 and delete them,
~ create new .csv files named 'MM_shLUC_{file_name}_shared_only.csv' and 'MM_WT_{file_name}_shared_only.csv',
~ files are created in 'Arkusze/{file_name} (Ekspr+MM)/Macierze Metylacji/Results',
(this files contain only shared genes between MM_df with localization in TSS1500 and TSS200,
------------------------

*s3MM_merge_files*
~ take 'MM_shLUC_{file_name}_shared_only.csv' and 'MM_WT_{file_name}_shared_only.csv',
~ create new .csv file with separate columns for probeID, localization, methylation value from WT and shLUC df,
~ created .csv file is named 'MM_{file_name}_shLUC_WT' in 'Arkusze/{file_name} (Ekspr+MM)/Macierze Metylacji/Results',
------------------------

*s4MM_DEG_merge*
~ take 'DEG_{file_name}' and 'MM_{file_name}_shLUC_WT',
~ find genes shared between this two files,
~ create new .csv file named 'Merg_{file_name}.csv' in 'Arkusze/{file_name} (Ekspr+MM)/Results',
------------------------

*s5Get_info*
~ take 'DEG_{file_name}' file, find: number of genes, number of up/down/down_up genes,
~ take 'MM_{file_name}_shLUC_WT' file, find: number of genes, number of hypermethylated/hypomethylated/hyper_hypomethylated genes,
~ take 'Merg_{file_name}' file, find: number of genes, number of: hyper_hypo|down, hyper_hypo|up, hyper|down, hyper|up, hypo|down, hypo|up,
------------------------
















