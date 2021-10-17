import numpy as np
import copy
import time
import os
import math
import sys
from scipy.optimize import minimize
from scipy.optimize import shgo

def ufiber_error(x):
    Em_bar = 0.506*1e4/1e6
    Em_inter = x[0]*1e4/1e6
    alpha = -x[1]
    num = x[2]

#    Em_bar = 5.06e3/1e6
#    Em_inter = 7.5426e3/1e6
#    alpha = -0.23465
#    num = 0.34

##################################################

    f=open('./source/umat.f90','r+')
    flist=f.readlines()

    flist[49]='    Em_bar = '+str(Em_bar) + '\n'
    flist[50]='    Em_inter = '+str(Em_inter) + '\n'
    flist[51]='    alpha = '+str(alpha) + '\n'
    flist[52]='    num = '+str(num) + '\n'

    f=open('./source/umat.f90','w+')
    f.writelines(flist)
    f.close()

    os.system('./compile.sh')

    os.system(solver_directory+' '+ jobname)
#    fiber_centroid_number = 68


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

    ufiber_error_disp = 0.0

    string = './iter.dat'
    f = open(string)
    iteration_step = int(f.readlines()[0])
    f.close()

    for i in np.arange(fiber_centroid_number):
        for meas_index in np.arange(measind1,measind2):
            ufiber_error_disp += ((ux[i] - uread[meas_index,2*i])**2 + (uy[i] - uread[meas_index,2*i+1])**2)**0.5/(uread[meas_index,2*i]**2 + uread[meas_index,2*i+1]**2)**0.5

    ufiber_error_force = abs((rf_ref-rf)/rf_ref)
    ufiber_error_disp = ufiber_error_disp/fiber_centroid_number/(measind2 - measind1)
    ufiber_error = np.log(ufiber_error_disp) + ufiber_error_force

    string = objdirectory_save + '/objfunc_e'+str(err_index)+\
            'Ei'+str(Eint_ind)+'a'+str(alpha_ind)+'num'+str(num_ind)+'.dat'
    f = open(string,'at')
    f.write(str(iteration_step)+' '+str(x[0])+' '+str(x[1])+' '+str(x[2])+'\n')
    f.write(str(iteration_step)+' '+str(ufiber_error_disp) + ' ' + str(ufiber_error_force) + ' ' + str(ufiber_error)+'\n')
    f.close()

    iteration_step += 1
    string = './iter.dat'
    f = open(string,'w')
    f.write(str(iteration_step))
    f.close()
    if iteration_step > 500:
        string = 'e'+str(err_index)+\
        'Ei'+str(Eint_ind)+'a'+str(alpha_ind)+'num'+str(num_ind) + ' Exceed maximum iteration steps(500)'
        sys.exit(string)
    return ufiber_error

#exact value [0.506, 0.75426, 0.23465, 0.34]

jobname = 'microns100'
solver_name = 'CalculiX'
err_index = 2
measind1 = 0
measind2 = 100
level = 4
Eint_ind = 2
alpha_ind = 2
num_ind = 2
cpu_ind = 1
fiber_centroid_number = 68
solver_directory1 = '/home/zsu/calculix/'
solver_directory2 = '/ccx_2.17/src/ccx_2.17'
solver_directory = solver_directory1 + solver_name + solver_directory2

parent_dir = '/home/zsu/NASAIRAD/opt_stratafiedsampling_estimationobjfunc/'
objdirectory_save = parent_dir + jobname + '_results/opt'
objdirectory_open = parent_dir + jobname + '_results/ref'

string = objdirectory_open+'/reactionforce.dat'
f = open(string)
lines = f.readlines()
rf_ref = float(lines[0])

for meas_index in np.arange(measind1,measind2):
    string = objdirectory_open+'/fiber_centroid_disp_e'+str(err_index)+'_m'+str(meas_index)+'.dat'
    f = open(string)
    lines = f.readlines()
    rows=len(lines)
    uread_temp = np.zeros((rows,2*fiber_centroid_number))
    row = 0
    for line in lines:
        list = line.strip('\n').split(' ')
        uread_temp[row,:] = list[::]
        row+=1
    if meas_index == measind1:
        uread = uread_temp
    else:
        uread = np.vstack([uread,uread_temp])
    f.close()

Einter_interval = np.arange(0.0,1.0 + 1.0/level,1.0/level)
alpha_interval = np.arange(0.0,1.0 + 1.0/level,1.0/level)
num_interval = np.arange(0.0,0.5 + 0.5/level,0.5/level)

Einter_interval[0] = Einter_interval[0] + 1e-2
alpha_interval[0] = alpha_interval[0] + 1e-2
num_interval[0] = num_interval[0] + 1e-2

Einter_interval[-1] = Einter_interval[-1] - 1e-2
alpha_interval[-1] = alpha_interval[-1] - 1e-2
num_interval[-1] = num_interval[-1] - 1e-2

bnds = [(0.0,1.0),(0.0,1.0),(0.0,0.499)]

x0_Einter = (Einter_interval[Eint_ind+1] - Einter_interval[Eint_ind]) * np.random.rand() + Einter_interval[Eint_ind]
x0_alpha = (alpha_interval[alpha_ind+1] - alpha_interval[alpha_ind]) * np.random.rand() + alpha_interval[alpha_ind]
x0_num = (num_interval[int(num_ind)+1] - num_interval[int(num_ind)]) * np.random.rand() + num_interval[int(num_ind)]
#x0_Einter = (Einter_interval[Eint_ind+1] - Einter_interval[Eint_ind]) * 0.5 + Einter_interval[Eint_ind]
#x0_alpha = (alpha_interval[alpha_ind+1] - alpha_interval[alpha_ind]) * 0.5 + alpha_interval[alpha_ind]
#x0_num = (num_interval[int(num_ind)+1] - num_interval[int(num_ind)]) * 0.5 + num_interval[int(num_ind)]
x0 = [x0_Einter, x0_alpha, x0_num]

string = objdirectory_save + '/objfunc_e'+str(err_index)+\
        'Ei'+str(Eint_ind)+'a'+str(alpha_ind)+'num'+str(int(num_ind))+'.dat'
f = open(string,'w')
f.close()

string = './iter.dat'
f = open(string,'w')
f.write(str(0))
f.close()

res = minimize(ufiber_error, x0, method='SLSQP', bounds=bnds, tol = 1e-6,\
               options={'maxiter': 1000,'disp': True,'ftol': 1e-6,'eps':1e-4})

print(res.x)
print(res.message)
print(res.fun)
