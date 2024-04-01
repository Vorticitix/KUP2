#import necessary packages
from glob import glob
import os
import sys

#set path and filenames 
path='/storage/climatestor/PleioCEP/doensen/data/'
files_30 = sorted(glob(path+'cyclone_era5_????_????/fort_30_????_????'))
files_36 = sorted(glob(path+'cyclone_era5_????_????/fort_36_????_????'))
#if output file already exists remove it
#try:
#    os.remove(path+"fort_30_total_RCP85.txt")
#except:
#    pass

try:
    os.remove(path+"fort_36_total_ERA5.txt")
except:
    pass

print(files_36)
#Append all fort.30 files into one txt file
#with open(path+"fort_30_total.txt",'a+') as foo:
#    for i,file in enumerate(files_30):
#        with open(file,'r') as read:
#            lines = read.readlines()
#            if i==0:
#                for line in lines:
#                    foo.write(line)
#            else:
#                for line in lines[2:]:
#                    foo.write(line)
#Append all fort.34 files into one txt file
with open(path+"fort_36_total_ERA5.txt",'a+') as foo:
    for i,file in enumerate(files_36):
        with open(file,'r') as read:
            lines = read.readlines()
            if i==0:
                for line in lines:
                    foo.write(line)
            else:
                for line in lines[2:]:
                    foo.write(line)  
    



