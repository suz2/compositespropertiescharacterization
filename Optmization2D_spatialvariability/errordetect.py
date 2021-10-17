import numpy as np
import os
import re
error_bnd = [0.0,1.0,2.0,4.0,8.0]
i = 0
for e_ind in np.arange(len(error_bnd)):
    if error_bnd[e_ind] == 0.0:
        cpu_ind = 1
        path = 'e'+ str(int(error_bnd[e_ind])) + '_cpu' + str(cpu_ind)
        filename_re = 'e'+ str(int(error_bnd[e_ind])) + '_cpu' + str(cpu_ind)
        for filename in os.listdir(path):
            if re.match(filename_re + ".e\d+", filename):
               f=open(path+'/'+filename,'r+')
               flist=f.readlines()
               # errorinfo = flist[-1].split()
               if len(flist) > 0:
                  for line in flist:
                     print(line.split())
    else:
        for cpu_ind in np.arange(1,5):
           path = 'e'+ str(int(error_bnd[e_ind])) + '_cpu' + str(cpu_ind)
           filename_re = 'e'+ str(int(error_bnd[e_ind])) + '_cpu' + str(cpu_ind)
           for filename in os.listdir(path):
              if re.match(filename_re + ".e\d+", filename):
                 f=open(path+'/'+filename,'r+')
                 flist=f.readlines()
                 # errorinfo = flist[-1].split()
                 if len(flist) > 0:
                    for line in flist:
                        print(line.split()[0])
