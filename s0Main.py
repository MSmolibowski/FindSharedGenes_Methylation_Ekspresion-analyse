#-----------------Description--------------------
#-----------------Libraries----------------------
import  subprocess
#----------------Functions-----------------------
#********************************************************************
#********************************************************************

#---------------CODE-----------------------------
file_name = input("Write file name ex: 'SKMES_sh714' (don't write .csv): ")         #get input from user

file = open(f'file_name.txt', "w")
file.write(file_name)                                                              #write file name to txt file
file.close()

print(f'START SCRIPTS FOR: {file_name} files.')
print("*******************************")
print("START: s1DEG_find_genes & s2MM_filter")
subprocess.run("python3 s1DEG_find_genes.py & python3 s2MM_filter.py", shell = True)    #run script 1 and 2 in the same time
print("END: s0Main")
print("*******************************")
print("START: s3MM_merge_files")
subprocess.run("python3 s3MM_merge_files.py", shell = True)                             #run script 3
print("END: s3MM_merge_files")
print("*******************************")
print("START: s4MM_DEG_merge")
subprocess.run("python3 s4MM_DEG_merge.py", shell = True)
print("END: s4MM_DEG_merge")
print("*******************************")
print("START: s5Get_info")
subprocess.run("python3 s5Get_info.py", shell = True)
print("END: s5Get_info")
print("*******************************")
print("END OF ANALYZE")