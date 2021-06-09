# Susceptible cost 1E5 with same death costs
# See where the line moves to


from simulation import *
from tbPlot import *
from datetime import datetime
import os.path 
import numpy as np

from initial import initial
import numpy as np 

if __name__ == '__main__':

    ########################################################
    ### TESTBENCH
    ########################################################

    dictVar = initial()

    nt = dictVar['nt']
    cVals = dictVar['cVals']

    #########
    ### RUN 1 
    #########

    # nsPerInt = [nt//4,nt//4,nt//4,nt- 3* nt//4]# no. of integration steps per interval


    # nIntervals = len(nsPerInt)
    # # Enumerative: discrete values on discrete intervals
    # ctrlVals = [ [ cVals[0,:], [], cVals[2,:], [] ]  for ii in range(nIntervals) ]
    # # Which controls are equal to each other
    # ctrlEq = [ [0,0,0,2]  for ii in range(nIntervals) ]     

    

    # deathArr = [ [ x, 0.75*x ] for x in 10 ** np.arange(4,6.7,0.1) ]
    # susceptibleArr =  [ [ x , x ] for x in 10 ** np.arange(0,6.7,0.2)]
    # susceptible_Death_Simulation(dictVar,nsPerInt,ctrlVals,ctrlEq,
    #                             deathArr,susceptibleArr,nProcesses = 10, filename = 'susceptible_death.json')   

    #########
    ### RUN 2 without modifying cost
    #########

    # check to see how good MC is

    # nsPerInt = [nt] #  no. of integration steps per interval


    # nIntervals = len(nsPerInt)
    # # Enumerative: discrete values on discrete intervals
    # ctrlVals = [ [ cVals[0,:], [], cVals[2,:], [] ]  for ii in range(nIntervals) ]
    # # Which controls are equal to each other
    # ctrlEq = [ [0,0,0,2]  for ii in range(nIntervals) ]     

    # deathArr = [ [ x, 0.75*x ] for x in 10 ** np.arange(4,6.7,0.5) ]
    # susceptibleArr =  [ [ x , x ] for x in 10 ** np.arange(0,6.7,0.5)]
    # susceptible_Death_Simulation(dictVar,nsPerInt,ctrlVals,ctrlEq,
    #                             deathArr,susceptibleArr,nProcesses = 10, filename = 'susceptible_death2.json')   

    #########
    ### RUN 3 with modifying cost
    #########

    # nsPerInt = [nt] #  no. of integration steps per interval

    # nIntervals = len(nsPerInt)
    # # Enumerative: discrete values on discrete intervals
    # ctrlVals = [ [ cVals[0,:], [], cVals[2,:], [] ]  for ii in range(nIntervals) ]
    # # Which controls are equal to each other
    # ctrlEq = [ [0,0,0,2]  for ii in range(nIntervals) ]     

    # deathArr = [ [ x, 0.75*x ] for x in 10 ** np.arange(4,6.6,0.2) ]
    # susceptibleArr =  [ [ x , x ] for x in 10 ** np.arange(0,6.6,0.5)]
    # susceptible_Death_Simulation(dictVar,nsPerInt,ctrlVals,ctrlEq,
    #                             deathArr,susceptibleArr,nProcesses = 10, filename = 'susceptible_death3.json')   


    # #########
    # ### RUN 4 ( NoImproveMax = 1 , flatIter = 10)
    # #########

    # nsPerInt = [nt] #  no. of integration steps per interval

    # nIntervals = len(nsPerInt)
    # # Enumerative: discrete values on discrete intervals
    # ctrlVals = [ [ cVals[0,:], [], cVals[2,:], [] ]  for ii in range(nIntervals) ]
    # # Which controls are equal to each other
    # ctrlEq = [ [0,0,0,2]  for ii in range(nIntervals) ]     

    # deathArr = [ [ x, 0.75*x ] for x in 10 ** np.arange(4.0,6.6,0.1) ]
    # susceptibleArr =  [ [ x , x ] for x in 10 ** np.arange(0.0,6.6,0.2)]
    # susceptible_Death_Simulation(dictVar,nsPerInt,ctrlVals,ctrlEq,
    #                             deathArr,susceptibleArr,nProcesses = 10, filename = 'susceptible_death4.json')   


    # #########
    # ### RUN 5 ( NO MC)
    # #########

    # nsPerInt = [nt] #  no. of integration steps per interval

    # nIntervals = len(nsPerInt)
    # # Enumerative: discrete values on discrete intervals
    # ctrlVals = [ [ cVals[0,:], [], cVals[2,:], [] ]  for ii in range(nIntervals) ]
    # # Which controls are equal to each other
    # ctrlEq = [ [0,0,0,2]  for ii in range(nIntervals) ]     

    # deathArr = [ [ x, 0.75*x ] for x in 10 ** np.arange(5.0,6.6,0.1) ]
    # susceptibleArr =  [ [ x , x ] for x in 10 ** np.arange(3.0,6.6,0.2)]
    # susceptible_Death_Simulation(dictVar,nsPerInt,ctrlVals,ctrlEq,
    #                             deathArr,susceptibleArr,nProcesses = 10, filename = 'susceptible_death5.json')  


    #########
    ### RUN 6 ( NO MC and uTol = 0.005 and other tolerance in ln578 set to 0.0005)
    #########

    # nsPerInt = [nt] #  no. of integration steps per interval

    # nIntervals = len(nsPerInt)
    # # Enumerative: discrete values on discrete intervals
    # ctrlVals = [ [ cVals[0,:], [], cVals[2,:], [] ]  for ii in range(nIntervals) ]
    # # Which controls are equal to each other
    # ctrlEq = [ [0,0,0,2]  for ii in range(nIntervals) ]     

    # deathArr = [ [ x, 0.75*x ] for x in 10 ** np.arange(5.0,6.6,0.2) ]
    # susceptibleArr =  [ [ x , x ] for x in 10 ** np.arange(3.0,6.6,0.5)]
    # susceptible_Death_Simulation(dictVar,nsPerInt,ctrlVals,ctrlEq,
    #                             deathArr,susceptibleArr,nProcesses = 5, filename = 'susceptible_death6.json')  


    # #########
    # ### RUN 7 ( NO ENUM NO MC and uTol = 0.005 and other tolerance in ln578 set to 0.0005)
    # #########

    # nsPerInt = [nt] #  no. of integration steps per interval

    # nIntervals = len(nsPerInt)
    # # Enumerative: discrete values on discrete intervals
    # ctrlVals = [ [ cVals[0,:], [], cVals[2,:], [] ]  for ii in range(nIntervals) ]
    # # Which controls are equal to each other
    # ctrlEq = [ [0,0,0,2]  for ii in range(nIntervals) ]     

    # deathArr = [ [ x, 0.75*x ] for x in 10 ** np.arange(5.0,6.6,0.2) ]
    # susceptibleArr =  [ [ x , x ] for x in 10 ** np.arange(3.0,6.6,0.5)]
    # susceptible_Death_Simulation(dictVar,nsPerInt,ctrlVals,ctrlEq,
    #                             deathArr,susceptibleArr,nProcesses = 5, filename = 'susceptible_death7.json')  



    #########
    ### RUN 8 ( Small runs for MC and uTol = 0.005 and other tolerance in ln578 set to 0.0005)
    #########

    # nsPerInt = [nt] #  no. of integration steps per interval

    # nIntervals = len(nsPerInt)
    # # Enumerative: discrete values on discrete intervals
    # ctrlVals = [ [ cVals[0,:], [], cVals[2,:], [] ]  for ii in range(nIntervals) ]
    # # Which controls are equal to each other
    # ctrlEq = [ [0,0,0,2]  for ii in range(nIntervals) ]     

    # deathArr = [ [ x, 0.75*x ] for x in 10 ** np.arange(5.0,6.6,0.2) ]
    # susceptibleArr =  [ [ x , x ] for x in 10 ** np.arange(3.0,6.6,0.5)]
    # susceptible_Death_Simulation(dictVar,nsPerInt,ctrlVals,ctrlEq,
    #                             deathArr,susceptibleArr,nProcesses = 10, filename = 'susceptible_death8.json')  



    ########
    ## nonImm_Death_partialEnum_partialMC:  partial enum , partial MC
    ########
    # print('---------------')
    # print('Imm_Death_partialEnum_partialMC')
    # print('---------------')

    # #  no. of integration steps per interval
    # nsPerInt = [nt//3, nt//3, nt - 2 * nt//3] 

    # nIntervals = len(nsPerInt)
    # # Enumerative: discrete values on discrete intervals
    # ctrlVals = [ [ cVals[0,:], [], cVals[2,:], [] ]  for ii in range(nIntervals) ]
    # # Which controls are equal to each other
    # ctrlEq = [ [0,0,0,2]  for ii in range(nIntervals) ]     

    # # Temperature Parameters
    # tempParam = {'flatIter':10,
    #              'noImproveMax':2}

    # deathArr = [ [ x, 0.75*x ] for x in 10 ** np.linspace(4.5, 6.5, 25) ]
    # susceptibleArr =  [ [ x , x ] for x in 10 ** np.linspace(0.0, 6.5, 25)]
    # susceptible_Death_Simulation(dictVar,nsPerInt,ctrlVals,ctrlEq,
    #                             deathArr,susceptibleArr,nProcesses = 10, 
    #                             filename = 'Imm_Death_partialEnum_partialMC.json',
    #                             tempParam = tempParam,
    #                             enumOn = True,
    #                             mcOn = True)  

    # ########
    # ## nonImm_Death_noEnum_fullMC:  No enum , Full MC
    # ########
    # print('---------------')
    # print('nonImm_Death_noEnum_fullMC')
    # print('---------------')

    # #  no. of integration steps per interval
    # nsPerInt = [] 
    # nIntervals = 0
    # # Enumerative: discrete values on discrete intervals
    # ctrlVals = []
    # # Which controls are equal to each other
    # ctrlEq = []     

    # # Temperature Parameters
    # tempParam = {'flatIter':10,
    #              'noImproveMax':5}

    # deathArr = [ [ x, 0.75*x ] for x in 10 ** np.linspace(4.5, 6.5, 25) ]
    # susceptibleArr =  [ [ x , x ] for x in 10 ** np.linspace(0.0, 6.5, 25)]
    # susceptible_Death_Simulation(dictVar,nsPerInt,ctrlVals,ctrlEq,
    #                             deathArr,susceptibleArr,nProcesses = 10, 
    #                             filename = 'nonImm_Death_noEnum_fullMC.json',
    #                             tempParam = tempParam,
    #                             enumOn = False,
    #                             mcOn = True)  


    ########
    ## nonImm_Death_fullEnum_noMC:  full enum , no MC
    ########
    print('---------------')
    print('nonImm_Death_fullEnum_noMC')
    print('---------------')

    #  no. of integration steps per interval
    nsPerInt = [nt//4, nt//4, nt//4, nt - 3 * nt//4]

    nIntervals = len(nsPerInt)
    # Enumerative: discrete values on discrete intervals
    ctrlVals = [ [ cVals[0,:], [], cVals[2,:], [] ]  for ii in range(nIntervals) ]
    # Which controls are equal to each other
    ctrlEq = [ [0,0,0,2]  for ii in range(nIntervals) ]     

    # Temperature Parameters
    tempParam = {'flatIter':10,
                 'noImproveMax':0}

    deathArr = [ [ x, 0.75*x ] for x in 10 ** np.linspace(4.5, 6.5, 25) ]
    susceptibleArr =  [ [ x , x ] for x in 10 ** np.linspace(0.0, 6.5, 25)]
    susceptible_Death_Simulation(dictVar,nsPerInt,ctrlVals,ctrlEq,
                                deathArr,susceptibleArr,nProcesses = 10, 
                                filename = 'nonImm_Death_fullEnum_noMC.json',
                                tempParam = tempParam,
                                enumOn = True,
                                mcOn = False)  
