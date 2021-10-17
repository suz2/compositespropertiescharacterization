import os
import sys
import numpy as np

def optdirectory(directory,job_name,err_index,measind1,measind2,solver_no,level,num_cpu,cpu_ind):

    os.system('mkdir '+directory)
    for Eint_ind in np.arange(level):
        for alpha_ind in np.arange(level):
            for num_ind in np.arange((cpu_ind - 1) * level / num_cpu, cpu_ind * level / num_cpu, 1):
                os.system('cp optmechload_expmod.py '+directory+'/optmechload_expmod_e'+str(err_index)+\
                            'Ei'+str(Eint_ind)+'a'+str(alpha_ind)+'num'+str(int(num_ind))+'.py')
    os.system('cp compile.sh '+directory)
    os.system('cp '+job_name+'.inp '+directory)
    os.system('cp UMATSRC.inc '+directory)
    os.system('cp -r source '+directory)
    os.system('cp initial.sh '+directory)
    ##############################
    f=open(directory+'/compile.sh','r+')
    flist=f.readlines()

    flist[13]='SOLVER_NAME=\'CalculiX_'+str(solver_no)+'\'\n'

    f=open(directory+'/compile.sh','w+')
    f.writelines(flist)
    f.close()
    ##############################
    for Eint_ind in np.arange(level):
        for alpha_ind in np.arange(level):
            for num_ind in np.arange((cpu_ind - 1) * level / num_cpu, cpu_ind * level / num_cpu, 1):
                f=open(directory+'/optmechload_expmod_e'+str(err_index)+\
                        'Ei'+str(Eint_ind)+'a'+str(alpha_ind)+'num'+str(int(num_ind))+'.py','r+')
                flist=f.readlines()
                flist[106]= 'jobname = \''+job_name+'\'\n'
                flist[107]= 'solver_name = \'CalculiX_'+str(solver_no)+'\'\n'
                flist[108]= 'err_index = '+str(err_index) + '\n'
                flist[109]= 'measind1 = '+str(measind1) + '\n'
                flist[110]= 'measind2 = '+str(measind2) + '\n'
                flist[111]= 'level = '+str(level) + '\n'
                flist[112]= 'Eint_ind = '+str(Eint_ind) + '\n'
                flist[113]= 'alpha_ind = '+str(alpha_ind) + '\n'
                flist[114]= 'num_ind = '+str(int(num_ind)) + '\n'
                flist[115]= 'cpu_ind = '+str(int(cpu_ind)) + '\n'

                f=open(directory+'/optmechload_expmod_e'+str(err_index)+\
                        'Ei'+str(Eint_ind)+'a'+str(alpha_ind)+'num'+str(int(num_ind))+'.py','w+')
                f.writelines(flist)
                f.close()
    ##############################
    f=open(directory+'/initial.sh','r+')
    flist=f.readlines()
    flist[1]= '#$ -N ' + 'e'+str(err_index) +'_cpu'+str(cpu_ind) + '\n'
    ii = 0
    time_line = '''duration=$SECONDS\necho "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."'''
    for Eint_ind in np.arange(level):
        for alpha_ind in np.arange(level):
            for num_ind in np.arange((cpu_ind - 1) * level / num_cpu, cpu_ind * level / num_cpu, 1):
                flist.append( '$HOME/python/Python-3.6.8/python ' + 'optmechload_expmod_e'+str(err_index)+\
                        'Ei'+str(Eint_ind)+'a'+str(alpha_ind)+'num'+str(int(num_ind))+'.py\n')
    flist.append(time_line)
    f=open(directory+'/initial.sh','w+')
    f.writelines(flist)
    f.close()
    return

job_name = 'microns100'
#solver_directory1 = '/home/zimu/MUMSresearch/calculix/'
solver_directory1 = '/home/zsu/calculix/'
#parent_dir = '/home/zimu/MUMSresearch/IRAD/opt_stratafiedsampling_clustertest/'
parent_dir = '/home/zsu/NASAIRAD/opt_stratafiedsampling_estimationobjfunc/'
#errorbnd = [0.0,0.01,0.02,0.04,0.08]
errorbnd = [0.0,0.01,0.02,0.04,0.08]
measind1 = 0
measind2 = 1
#stratafied_level
level = 4
# num_cpu should be even number
num_cpu = 4
fiber_centroid_number = 197

os.system('mkdir '+job_name+'_results')
os.system('mkdir '+job_name+'_results/ref')
os.system('mkdir '+job_name+'_results/opt')

os.system('rm -r e*_cpu*')

##############################
f=open('optmechload_expmod.py','r+')
flist=f.readlines()
flist[116]= 'fiber_centroid_number = ' + str(fiber_centroid_number) + '\n'
flist[117]= 'solver_directory1 = \'' + solver_directory1 + '\'\n'
flist[121]= 'parent_dir = \'' + parent_dir + '\'\n'
f=open('optmechload_expmod.py','w+')
f.writelines(flist)
f.close()

##############################
f=open('compile.sh','r+')
flist=f.readlines()
flist[14]= 'SOLVER_DIR1=' + '\'' + solver_directory1 + '\'' + '\n'
f=open('compile.sh','w+')
f.writelines(flist)
f.close()

##############################
f=open('source/umat.f90','r+')
flist=f.readlines()
flist[55]= '    dir = \'' + parent_dir + '\'\n'
f=open('source/umat.f90','w+')
f.writelines(flist)
f.close()

solver_no = 1
f_sol = open('solver_no_info.dat','w')

for index in np.arange(len(errorbnd)):
    err_index = int(errorbnd[index]*100)
    if err_index == 0:
        for cpu_ind in (np.arange(num_cpu) + 1):
            directory = 'e'+str(err_index)+'_cpu'+str(cpu_ind)
            optdirectory(directory,job_name,err_index,0,1,solver_no,level,num_cpu,cpu_ind)
            f_sol.write('solver_no '+str(solver_no)+' ' + directory + '\n')
            solver_no += 1
            ##os.system('mkdir '+job_name+'_results/opt/'+'e'+str(err_index)+'_cpu'+str(cpu_ind))
    else:
        for cpu_ind in (np.arange(num_cpu) + 1):
            directory = 'e'+str(err_index)+'_cpu'+str(cpu_ind)
            optdirectory(directory,job_name,err_index,measind1,measind2,solver_no,level,num_cpu,cpu_ind)
            f_sol.write('solver_no '+str(solver_no)+' ' + directory + '\n')
            solver_no += 1
            ##os.system('mkdir '+job_name+'_results/opt/'+'e'+str(err_index)+'_cpu'+str(cpu_ind))

f_sol.close()
