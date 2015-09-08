############################################
###
#
#   This file contains codes for processing
#   data of density of states, usually got
#   from Wang-Landau simulation
#
#   @Author: Jerry Shi
#   @Date: Jun 6th, 2013
#
###
############################################


import numpy as np

#####################
# @input:
#   1. ndarray: [0] energy [1-] ln(dos) for runs
#   2. dt
#   3. num_points (temperature points)
#   4. sys_size (size of system, to calculate specific heat)
# @return:
#   1. ndarray: [0] temperature [1] Cv [2] eb_Cv [3] aveE [4] eb_aveE
def thermo(ln_dos, dT, num_points, sys_size):
    kb=1
    num_runs = ln_dos.shape[1] -1
    thermo_result=np.zeros((num_points,5))
    size = ln_dos.shape[0] # number of points of density of states

    for i in range(0,num_points):
        T = i*dT+dT
        sumE = np.zeros((1,num_runs))
        sumE2 = np.zeros((1,num_runs))
        factor = np.reshape(-ln_dos[0,0]/(kb*T),(1,1))
        normal = np.reshape(ln_dos[0,1:num_runs+1], (1,num_runs)) + factor
        for j in range(1,size):
            factor = np.reshape(-ln_dos[j,0]/(kb*T),(1,1))
            normal += np.reshape( np.log1p( np.exp( np.reshape(ln_dos[j,1:num_runs+1],(1,num_runs)) + factor - normal) )  ,(1,num_runs))

        for j in range(0,size):
            factor = np.reshape(-ln_dos[j,0]/(kb*T),(1,1))
            poss = np.reshape( np.exp( np.reshape(ln_dos[j,1:num_runs+1],(1,num_runs)) + factor - normal),(1,num_runs))
            sumE  += ln_dos[j,0]*poss
            sumE2 += np.square(ln_dos[j,0])*poss

        sumE = np.reshape(sumE, (1,num_runs))
        sumE2 = np.reshape(sumE2, (1,num_runs))
        cv = np.reshape((sumE2-sumE*sumE)/T/T/sys_size,(1,num_runs))

        thermo_result[i] =  [   T                               ,
                                np.mean(cv,axis=1)              ,
                                np.std(cv,axis=1)               ,
                                np.mean(sumE,axis=1)            ,
                                np.std(sumE,axis=1)             ]
    return thermo_result

