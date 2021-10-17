import numpy as np
import os
import re
path = "./stratafied_sampling_2level"

error_bnd = [0.0,1.0,2.0,4.0]

time = np.array([[0,0]])
i = 0
for e_ind in np.arange(len(error_bnd)):
    if error_bnd[e_ind] == 0:
        m_ind = 0
        filename_re = 'e'+ str(int(error_bnd[e_ind])) + 'm' + str(m_ind)

        for filename in os.listdir(path):
            if re.match(filename_re + ".o\d+", filename):
                f=open(path+'/'+filename,'r+')
                flist=f.readlines()
                timeline = flist[-1].split()
                if i == 0:
                    time = np.array([[error_bnd[e_ind],float(timeline[0]) + float(timeline[3])/60.0]])
                else:
                    time = np.concatenate((time,np.array([[error_bnd[e_ind],float(timeline[0]) + float(timeline[3])/60.0]])))

                # time[i,0] = error_bnd[e_ind]
                # time[i,1] = float(timeline[0]) + float(timeline[3])/60.0
                i += 1
    else:
        for m_ind in np.arange(5):
            filename_re = 'e'+ str(int(error_bnd[e_ind])) + 'm' + str(m_ind)

            for filename in os.listdir(path):
                if re.match(filename_re + ".o\d+", filename):
                    f=open(path+'/'+filename,'r+')
                    flist=f.readlines()
                    timeline = flist[-1].split()
                    if i == 0:
                        time = np.array([[error_bnd[e_ind],float(timeline[0]) + float(timeline[3])/60.0]])
                    else:
                        time = np.concatenate((time,np.array([[error_bnd[e_ind],float(timeline[0]) + float(timeline[3])/60.0]])))
                    i += 1
print(time)
print(len(time))
print(time[0])
f = open(path + '/timecomputing.dat','w')
for i in np.arange(len(time)):
    f.write(str(time[i][0]) + ' ' + str(time[i][1]) + '\n')

f.close()