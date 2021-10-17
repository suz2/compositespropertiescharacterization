import numpy as np
import copy
import time
import sys,os
import math

os.system('./compile.sh')
jobname = 'microns200vf42'
solver_name = 'CalculiX'
solver_directory1 = '/home/zimu/MUMSresearch/calculix/'
solver_directory2 = '/ccx_2.17/src/ccx_2.17'
solver_directory = solver_directory1 + solver_name + solver_directory2

parent_dir = '/home/zimu/MUMSresearch/IRAD/opt_stratafiedsampling_clustertest/'
directory = parent_dir + jobname + '_results/ref'

measind1 = 0
measind2 = 100

errorbnd = [0.0,0.01,0.02,0.04,0.08]

Em_bar = 5.06e3/1e6
Em_inter = 7.5426e3/1e6
alpha = -0.23465
num = 0.34

f=open('./source/umat.f90','r+')
flist=f.readlines()

flist[49]='    Em_bar = '+str(Em_bar) + '\n'
flist[50]='    Em_inter = '+str(Em_inter) + '\n'
flist[51]='    alpha = '+str(alpha) + '\n'
flist[52]='    num = '+str(num) + '\n'

f=open('./source/umat.f90','w+')
f.writelines(flist)
f.close()

os.system(solver_directory+' '+ jobname)
fiber_centroid_number = 197

f = open('./' + jobname + '.dat', 'r')
lines = f.readlines()
displines = lines[3:3+fiber_centroid_number]

node_ind = 0
ux = np.zeros(fiber_centroid_number)
uy = np.zeros(fiber_centroid_number)
for line in displines:
    displist = line.strip('\n').split(' ')
    disp = []
    for t in displist:
        try:
            disp.append(float(t))
        except ValueError:
            pass
    ux[node_ind] = disp[1]
    uy[node_ind] = disp[2]
    node_ind += 1

print('ux:',ux)
print('uy:',uy)

rflines = lines[6+fiber_centroid_number:]
rf = 0
for line in rflines:
    rflist = line.strip('\n').split(' ')
    rflist2 = []
    for t in rflist:
        try:
            rflist2.append(float(t))
        except ValueError:
            pass
    rf += rflist2[2]
print('rf:',rf)

f = open(directory+'/reactionforce.dat','w')
f.write(str(rf))
f.close()

for index in np.arange(len(errorbnd)):
    err_index = int(errorbnd[index]*100)

    if err_index == 0:
        meas_ind1 = 0
        meas_ind2 = 1
    else:
        meas_ind1 = measind1
        meas_ind2 = measind2

    for ii in np.arange(meas_ind1,meas_ind2,1):
        f = open(directory+'/fiber_centroid_disp_e'+str(err_index)+'_m'+str(ii)+'.dat','w')
        f_err_ux = open(directory+'/error_ux_e'+str(err_index)+'_m'+str(ii)+'.dat','w')
        f_err_uy = open(directory+'/error_uy_e'+str(err_index)+'_m'+str(ii)+'.dat','w')
        for j in np.arange(fiber_centroid_number):
             if j == fiber_centroid_number-1:
                error_ux = np.random.rand()*errorbnd[index]*2-errorbnd[index]
                error_uy = np.random.rand()*errorbnd[index]*2-errorbnd[index]
                f.write(str(ux[j]*(1+error_ux))+' ')
                f.write(str(uy[j]*(1+error_uy)))
                f_err_ux.write(str(error_ux))
                f_err_uy.write(str(error_uy))
             else:
                error_ux = np.random.rand()*errorbnd[index]*2-errorbnd[index]
                error_uy = np.random.rand()*errorbnd[index]*2-errorbnd[index]
                f.write(str(ux[j]*(1+error_ux))+' ')
                f.write(str(uy[j]*(1+error_uy))+' ')
                f_err_ux.write(str(error_ux)+' ')
                f_err_uy.write(str(error_uy)+' ')
        f.close()
        f_err_ux.close()
        f_err_uy.close()

        f = open(directory+'/reactionforce_e'+str(err_index)+'_m'+str(ii)+'.dat','w')
        f_err = open(directory+'/error_reactionforce_e'+str(err_index)+'_m'+str(ii)+'.dat','w')
        error = np.random.rand()*errorbnd[index]*2-errorbnd[index]
        f.write(str(rf*(1+error)))
        f_err.write(str(error))

        f.close()
        f_err.close()
