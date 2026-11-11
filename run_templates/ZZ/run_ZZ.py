import os
import sys
import shutil
#import numpy as np
import re
from glob import glob

def run():
    n_stat = 15 # number of random points is 2**(n_stat)

    i_seed = 1 # choose a random seed for MC sampling

    # for scale-uncertainty calculation,
    # choose factorization and renormalization scale factors,
    # 1 is the central-scale value:
    i_scale_F = 1.0 # factorization scale factor
    i_scale_R = 1.0 # renormalization scale factor

    isubmit = 0 # (0: run the executable directly on this machine; 1: submit to a cluster using sge_qsub;)

    # choose which calculation,
    # "off_shell", "POLAR_LO", "POLAR_NLOQCD", "POLAR_NLOQCDEW":
    CALCULATION = "POLAR_NLOQCDEW"

    # sub-dominant contributions for POLAR calculations:
    gg_process = True
    bb_process = 2 # (choose 0 for NO bb, 1 for LO, 2 for NLO)

    # re-compile amplitude dir (default: 0 [no re-compile], 1 [force re-compile], 2 [clean-up]):
    recompile_proc = 0

    # choose an output diretory name:
    postfix=''
    output_dir='results_'+str(n_stat)+postfix
    
    energy = 13000 # collision energy in GeV
    
    IDp_range = fill_IDp_range(CALCULATION,gg_process,bb_process)
    print("IDp_range =",IDp_range)

    IDc_range_all = fill_IDc_range(IDp_range,CALCULATION,gg_process,bb_process)
    print("IDc_range_all =",IDc_range_all)

    # polarizations:
    IDh_range = [22,55,52,25,44]
    #IDh_range = [22]
    if (CALCULATION == "off_shell"):
        IDh_range = [44] # only unpolarized cross section available

    # optimize scale factors (LO and NLOEW have only mu_F; NLOQCD both mu_F and mu_R)
       # re-calculate i_scale_F and i_scale_R depending on IDc
    # optimize statistics (more stats for large contributions)
       # re-calculate n_stat depending on IDc
    # more input parameters in the directory cards/

# checking a specific process:
#    IDp_range = [8]

    i_process = -1
    for IDp in IDp_range:
        i_process += 1
        IDc_range = IDc_range_all[i_process]

# checking a specific contribution:
#        IDc_range = [1]

        for IDc in IDc_range:
            for IDh in IDh_range:
                submit(output_dir,recompile_proc,isubmit,IDp,IDc,IDh,i_seed,i_scale_F,i_scale_R,n_stat,energy)

##############################################################
### Location of MulBos:
home_mulbos='/hdd1/ldninh/hahn/working/MulBos_public/github/MulBos'

### Location of LHAPDF:
LHAPDF = "/home/ldninh/codes_hep/lib/lhapdf-6.5.6"

home=os.getcwd()

mulbos_version_file=home_mulbos+'/mulbos_version.txt'

lhapdf_lib_file=LHAPDF + '/lib/libLHAPDF.a'

if os.path.isfile(mulbos_version_file):
    print("MulBos path OK:",home_mulbos)
else:
    print("MulBos path NOT OK:",home_mulbos)
    exit()

if os.path.isfile(lhapdf_lib_file):
    print("LHAPDF path OK:",LHAPDF)
else:
    print("LHAPDF path NOT OK:",LHAPDF)
    exit()

# adding MulBos home to the system path
#sys.path.insert(0, home_mulbos)

# Add the directory containing the processes_list.py to the python-search path
sys.path.append(home_mulbos + "/processes")

import processes_list
proce = processes_list.process

#print("process full =",proce)
#for i in range(100):
#    print("process =",i,proce[i])

LoopInts_old = "LT = $(HOME1)/tools/loops/LoopInts/include"
LoopInts_new = LoopInts_old

LoopQ = "1" # 1 (LoopInts), 2 (LoopInts)

#print("home dir=",home)
#print("mulbos dir=",home_mulbos)

# to replace compiler:
Compiler_old = "FC = ifx"
Compiler_new = "FC = ifx"

def fill_IDp_range(CALCULATION,gg_process,bb_process):
    IDp_range = [8,9]
    if (CALCULATION == "off_shell"):
        IDp_range += [11,12] # include bb and aa by default
    else: # polarization calculations
        if (gg_process):
            IDp_range += [10]
        if (bb_process > 0):
            IDp_range += [11]
    return IDp_range

def fill_IDc_range(IDp_range,CALCULATION,gg_process,bb_process):
# processes: UU_DD_gg_bb_aa = [8,9,10,11,12]
    n_process = len(IDp_range)
    n_UU_DD = 2 # 8,9

    gg_row = [1] #10
    bb_row_LO = [1] #11
    bb_row_NLOQCD = list(range(1, 7 + 1))
    bb_row_NLOQCDEW = list(range(1, 17 + 1))
    if (bb_process==1): # 11
        bb_row = bb_row_LO
    elif (bb_process==2 and CALCULATION=="POLAR_NLOQCD"):
        bb_row = bb_row_NLOQCD
    elif (bb_process==2 and CALCULATION=="POLAR_NLOQCDEW"):
        bb_row = bb_row_NLOQCDEW
    
    if (CALCULATION == "off_shell"):
        A = [[0] for _ in range(n_process)]
    else: # polarization calculations
        if (CALCULATION == "POLAR_LO"):
            A = [[1] for _ in range(n_process)]
        if (CALCULATION == "POLAR_NLOQCD"):
            A = [list(range(1, 7 + 1)) for _ in range(n_UU_DD)] # 8,9
            if (gg_process): # 10
                A.append(gg_row)
            if (bb_process>0): # 11
                A.append(bb_row)
        if (CALCULATION == "POLAR_NLOQCDEW"):
            A = [list(range(1, 17 + 1)) for _ in range(n_UU_DD)]
            if (gg_process): # 10
                A.append(gg_row)
            if (bb_process>0): # 11
                A.append(bb_row)
    return A

def replace_text(file_in,file_out,old_text,new_text):
    os.system("sed" + " 's|"+old_text+"|"+new_text+"|g'" + " " + file_in + " > " + file_out)

def submit(output_dir,recompile_proc,isubmit,IDp,IDc,IDh,i_seed,i_scale_F,i_scale_R,n_stat,energy):
    i_proc = IDp
    i_comp = IDc

    if (IDc==0): # full off-shell
        IDh = 444

    polar=str(IDh) # polarization

    j_proc_text=str(IDp)
    i_comp_text=str(i_comp)
    postfix='_'+j_proc_text+'_'+i_comp_text+'_pol'+polar+'_seed'+str(i_seed)+'_scale'+str(i_scale_F)+'_'+str(i_scale_R)
    N_comps = len(proce[i_proc])

    #print("process =",i_proc,proce[i_proc])
    #print("N_comps =",N_comps)

    if i_comp in range(N_comps):
        proc=proce[i_proc][i_comp]
    else:
        print("invalid process id or contribution id; IDp, IDc =",IDp, IDc)
        return()

    print("process=",proc)

    process = home_mulbos+'/processes/'+proc

    if (recompile_proc == 2):
        os.chdir(process)
        os.system("make clean")
        os.chdir(home)
        return()

    process_lib_file = process + '/' + 'libsquaredme.a'
    if (not os.path.isfile(process_lib_file) or (recompile_proc == 1)):
        os.chdir(process)
        replace_text("./makefile","./makefile.tmp",Compiler_old,Compiler_new)
        if (LoopQ=="2"):
            replace_text("./makefile.tmp","./makefile.tmp",LoopInts_old,LoopInts_new)
        os.system("make -f makefile.tmp")
        os.system("rm -f makefile.tmp")
        os.chdir(home)

    n_stat_in = n_stat
        
    print("n_stat_in =",n_stat_in)

    output = output_dir + '/'    
    output_run=output+'run'+postfix
    output_dist1='dist'
    output_dist2=output_run+'/'+output_dist1
    output_print='xs'+'.txt'
    check_file = os.path.isfile(output_dist2+'/dist401')
    if (check_file):
        return()
    else:
        os.system("mkdir -p %s" % output_dist2)
###[edit input files
# makefile:
    with open(home_mulbos+'/driver/makefile.in', 'r') as file:
        data = file.read()
        data = data.replace("PROC_INPUT", proc)
        data = data.replace("HOME2_INPUT", home_mulbos)
        data = data.replace("LHAPDF_INPUT", LHAPDF)
        data = data.replace("LoopQ_INPUT", LoopQ)
        data = data.replace(Compiler_old, Compiler_new)
    with open(output_run+'/makefile', 'w') as file:
        file.write(data)
# input cards:
# input parameters:
    with open(r'cards/input.dat.in', 'r') as file:
        data = file.read()
        data = data.replace('PROC_ID_INPUT', j_proc_text)
        data = data.replace('PROC_CONT_INPUT', i_comp_text)
        data = data.replace('RES_INPUT', output_dist1)
        data = data.replace('ECM_INPUT', str(energy)+"d0")
        data = data.replace('STAT_INPUT', str(n_stat_in))
        data = data.replace('POLAR_INPUT', polar)
        data = data.replace('SEED_INPUT', str(i_seed))
        data = data.replace('MUF_FAC_INPUT', str(i_scale_F)+"d0")
        data = data.replace('MUR_FAC_INPUT', str(i_scale_R)+"d0")
    with open(output_run+'/input.dat', 'w') as file:
        file.write(data)
# input cuts:
    with open(r'cards/cuts.dat.in', 'r') as file:
        data = file.read()
    with open(output_run+'/cuts.dat', 'w') as file:
        file.write(data)
# copy main fortran file:
    os.system("cp %s" % home_mulbos+"/driver/main.F " + output_run)
###]edit input files

###[compile and submit
    os.chdir(output_run)
    os.system("make")

    if (isubmit==0):
        os.system("./mulbos %s" % " > "+output_print+"&")
    elif (isubmit==1):
        qsub_check = shutil.which("qsub")
        if qsub_check is None:
            raise RuntimeError("SGE qsub is not on this computer")
        os.system("qsub "+home_mulbos+"/driver/sge_submit.sh")
    else:
        print("invalid value for isubmit:",isubmit)
        os.chdir(home)
        exit()
    os.chdir(home)
###]compile and submit

### run MulBos:
if __name__ == "__main__":
    run()


