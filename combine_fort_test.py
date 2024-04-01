#import necessary packages
from glob import glob
import os
import sys

#set path and filenames 
path='/storage/climatestor/PleioCEP/doensen/data/'
#files_30 = sorted(glob(path+'cyclone_*/fort_30_00[0123]?_00[1234]?'))
files_36 = sorted(glob(path+'cyclone_*/fort_36_00[34]0_00[34]9'))
#if output file already exists remove it
#files_36 = [path+'cyclone_0005_0604/fort_36_0010_0019',path+'cyclone_0005_0604/fort_36_0020_0029']
try:
    #os.remove(path+"fort_30_total_test_med.txt")
    os.remove(path+"fort_36_total_test_med.txt")
except:
    pass
print(files_36)
#Append all fort.30 files into one txt file
# with open(path+"fort_30_total_test.txt",'a+') as foo:
#     for i,file in enumerate(files_30):
#         with open(file,'r') as read:
#             lines = read.readlines()
#             if i==0:
#                 for line in lines:
#                     foo.write(line)
#             else:
#                 for line in lines[2:]:
#                     foo.write(line)
#Append all fort.34 files into one txt file
with open(path+"fort_36_total_test_med.txt",'a+') as foo:
    for i,file in enumerate(files_36):
        with open(file,'r') as read:
            lines = read.readlines()
            if i==0:
                for line in lines:
                    foo.write(line)
            else:
                for line in lines[2:]:
                    foo.write(line)  



