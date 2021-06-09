import numpy as np
import time
import json
from multiprocessing import Pool

from COVIDcodeMC import optMinCost


def cost_multiplier(dictVar,nsPerInt,ctrlVals,ctrlEq,sim = '',multiplierArr = [1], filename = None):
    '''
    Return:
        dictionary with the information collected in the simulation.
    Parameter:
        dictVar: Variables defined in initial.py
        nsPerInt: An array contains the number of steps per interval
        ctrlVals:
        ctrlEq:
        sim: A single character string that determines which simulation will be used [a,b,c,d,e,f]
        multiplier: An array of values that is used to multiply the cost
        filename: a string for the name of a json file.
    
    '''
    
    a = dictVar['a']
    b = dictVar['b']
    c = dictVar['c']
    d = dictVar['d']
    e = dictVar['e']
    f = dictVar['f']
    g = dictVar['g']
    params = dictVar['params']
    tf = dictVar['tf']
    nt = dictVar['nt']
    xI = dictVar['xI']
    cVals = dictVar['cVals']
    uMax = dictVar['uMax'] 


    # number of runs 
    nRuns = len(multiplierArr)
    
    # Initialize the arrays:  costArray (size is 8, number of runs)
    costArray = np.zeros( (nRuns, 9) )

    costBestArray = np.zeros(nRuns)
    
    uCtrl = np.zeros( (nRuns,tf,4) )

    timeArr = np.zeros(nRuns)


    for ii,multiplier in enumerate(multiplierArr):
        print('{} Simulation: {} out of {}'.format(sim,ii+1,nRuns))
        t1 = time.perf_counter()


        if sim == 'a':
            a = dictVar['a'] * multiplier
            a[:,0] = dictVar['a'][:,0]
        elif sim =='b':
            b = dictVar['b'] * multiplier
        elif sim =='c':
            c[0] = dictVar['c'][0] * multiplier
            c[1] = dictVar['c'][1] * multiplier
        elif sim =='d':
            d[0] = dictVar['d'][0] * multiplier
            d[1] = dictVar['d'][1] * multiplier
        elif sim =='e':
            e[0] = multiplier # dictVar['e'][0] * multiplier
            e[1] = 0.75 * e[0]           
        elif sim =='f':
            f[0] = dictVar['f'][0] * multiplier
            f[1] = dictVar['f'][1] * multiplier
        else:
            raise('ERROR: Simulation not found.')


        _,uBest,costBest,allCostsSummed = optMinCost(a,b,c,d,e,f,g,params,tf,nt,xI,nsPerInt,cVals,ctrlVals,ctrlEq,uMax)
    
        uCtrl[ii,:,:] = uBest
        costArray[ii,:] = allCostsSummed 
        costBestArray[ii] = costBest
    
        t2 = time.perf_counter()
        timeArr[ii] = t2-t1
        
   
    x = {'multiplier':multiplierArr.tolist(),
        'costArray':costArray.tolist(),
        'runtime':timeArr.tolist(),
        'costBest':costBestArray.tolist(),
        'uCtrl':uCtrl.tolist()}  
 

    #######################################
    ### SAVE DATA TO A JSON FILE
    #######################################  
    
    if filename is not None:

        with open(filename, 'w') as f:
            json.dump(x, f,indent = 4)

    return x


def func_optMinCost(runID,lstVar):

    try:
        a,b,c,d,e,f,g,params,tf,nt,xI,nsPerInt,cVals,ctrlVals,ctrlEq,uMax,kwargs = lstVar        
        # xBest,uBest,costBest,allCostsSummed 
        print('RunID:',runID)
        return optMinCost(a,b,c,d,e,f,g,params,tf,nt,xI,nsPerInt,cVals,ctrlVals,ctrlEq,uMax,**kwargs)

    except:

        print('FAILED')


def multi_processing(func,arr,processes = 1):

    runIDs = list(range(len(arr)))

    with Pool(processes= processes) as p:
        results = p.starmap(func , zip(runIDs,arr) )

    return results


def susceptible_Death_Simulation(dictVar,nsPerInt,ctrlVals,ctrlEq,deathArr,susceptibleArr,nProcesses = 1,filename = None,**kwargs):

    a = dictVar['a']
    b = dictVar['b']
    c = dictVar['c']
    d = dictVar['d']
    #e = dictVar['e']
    f = dictVar['f']
    #g = dictVar['g']
    params = dictVar['params']
    tf = dictVar['tf']
    nt = dictVar['nt']
    xI = dictVar['xI']
    cVals = dictVar['cVals']
    uMax = dictVar['uMax'] 

    # number of runs 
    nDeathRuns = len(deathArr)
    nSusRuns = len(susceptibleArr)

    print('Total Number of Runs:',nDeathRuns * nSusRuns)

    arr = [ [] for ii in range( nDeathRuns * nSusRuns )  ]
    kk = 0
    for e in deathArr:
        for g in susceptibleArr:

            arr[kk] = [a,b,c,d,e,f,g,params,tf,nt,xI,nsPerInt,cVals,ctrlVals,ctrlEq,uMax,kwargs]
            kk += 1


    results = multi_processing(func_optMinCost, arr, processes= nProcesses)
    


    # TOTAL COST FOR EVERYTHING
    enum_costBestArray = np.zeros( (nDeathRuns,nSusRuns) )
    mc_costBestArray = np.zeros( (nDeathRuns,nSusRuns) )
    opt_costBestArray = np.zeros( (nDeathRuns,nSusRuns) )
    # CONTROL MEASURES
    enum_uCtrl = np.zeros( (nDeathRuns,nSusRuns,tf,4) )
    mc_uCtrl = np.zeros( (nDeathRuns,nSusRuns,tf,4) )
    opt_uCtrl = np.zeros( (nDeathRuns,nSusRuns,tf,4) )


    # COST FOR THE COMPARTMENTS
    costArray = np.zeros( (nDeathRuns, nSusRuns,8) )
    # POPULATION 
    xBestArray = np.zeros( (nDeathRuns,nSusRuns, tf + 1 , 18) )
    #RUNTIME
    timeArr = np.zeros( (nDeathRuns,nSusRuns,3) )

    # INFORMATION ON THE RUNS IN MONTE CARLOS   
    iterCost = np.zeros( (nDeathRuns,nSusRuns) ).tolist()
    minCost = np.zeros( (nDeathRuns,nSusRuns) ).tolist()

    idxE = 0
    idxG = 0

    for returns in results:

        enumArr , mcArr , optArr , iterInfo , runtimes = returns 

        #xBest,uBest,costBest,allCostsSummed, runtimes = returns

        if idxG == nSusRuns:

            idxE += 1
            idxG = 0


        enum_uCtrl[idxE,idxG,:,:] = enumArr[0]
        mc_uCtrl[idxE,idxG,:,:] = mcArr[0]
        opt_uCtrl[idxE,idxG,:,:] = optArr[0]

        enum_costBestArray[idxE,idxG] = enumArr[1]
        mc_costBestArray[idxE,idxG] = mcArr[1]
        opt_costBestArray[idxE,idxG] = optArr[1]

        
        xBestArray[idxE,idxG,:] = optArr[2]     
        costArray[idxE,idxG,:] = optArr[3]
        timeArr[idxE,idxG,:] = runtimes

        minCost[idxE][idxG] = iterInfo[0]
        iterCost[idxE][idxG] = iterInfo[1]

        idxG += 1


    x = {'Parameter':{'uMax':uMax.tolist(),
                      'tf':tf},
         'changes':{'eArr':deathArr,
                   'gArr':susceptibleArr},
         'uCtrl':{'enum':enum_uCtrl.tolist(),
                  'mc':mc_uCtrl.tolist(),
                  'opt':opt_uCtrl.tolist()},
         'totCost':{'enum':enum_costBestArray.tolist(),
                    'mc':mc_costBestArray.tolist(),
                    'opt':opt_costBestArray.tolist()},
         'MCInfo':{'minCost':minCost,
                   'iterCost':iterCost},
        'costArr':costArray.tolist(),
        'xIArr': xBestArray.tolist(),
        'runtime':timeArr.tolist()}


    # x = {'changes':{'deathArr':deathArr,
    #                 'susceptibleArr':susceptibleArr},
    #     'xBestArray': xBestArray.tolist(),
    #     'costArray':costArray.tolist(),
    #     'runtime':timeArr.tolist(),
    #     'costBest':costBestArray.tolist(),
    #     'uCtrl':uCtrl.tolist()}  
 

    #######################################
    ### SAVE DATA TO A JSON FILE
    #######################################  
    
    if filename is not None:

        with open(filename, 'w') as f:
            json.dump(x, f,indent = 4)

    return x















    # a = dictVar['a']
    # b = dictVar['b']
    # c = dictVar['c']
    # d = dictVar['d']
    # e = dictVar['e']
    # f = dictVar['f']
    # g = dictVar['g']
    # params = dictVar['params']
    # tf = dictVar['tf']
    # nt = dictVar['nt']
    # xI = dictVar['xI']
    # cVals = dictVar['cVals']
    # uMax = dictVar['uMax'] 



    # # number of runs 
    # nDeathRuns = len(deathArr)
    # nSusRuns = len(susceptibleArr)
    
    # # Initialize the arrays:  costArray (size is 9, number of runs)
    # costArray = np.zeros( (nDeathRuns, nSusRuns,9) )
    # costBestArray = np.zeros( (nDeathRuns,nSusRuns) )
    # uCtrl = np.zeros( (nDeathRuns,nSusRuns,tf,4) )
    # timeArr = np.zeros( (nDeathRuns,nSusRuns) )

    # nRuns = nDeathRuns * nSusRuns
    # kk = 1

    # for ii,e in enumerate(deathArr):
    #     for jj,g in enumerate(susceptibleArr):
    #         print('{} out of {} Runs'.format(kk,nRuns))

    #         t1 = time.perf_counter()

    #         _,uBest,costBest,allCostsSummed = optMinCost(a,b,c,d,e,f,g,params,tf,nt,xI,nsPerInt,cVals,ctrlVals,ctrlEq,uMax)
        
    #         uCtrl[ii,jj,:,:] = uBest
    #         costArray[ii,jj,:] = allCostsSummed 
    #         costBestArray[ii,jj] = costBest
        
    #         t2 = time.perf_counter()
    #         timeArr[ii,jj] = t2-t1      

    #         kk += 1   


    # x = {'deathArr':deathArr,
    #     'susceptibleArr':susceptibleArr,
    #     'costArray':costArray.tolist(),
    #     'runtime':timeArr.tolist(),
    #     'costBest':costBestArray.tolist(),
    #     'uCtrl':uCtrl.tolist()}  
 

    # #######################################
    # ### SAVE DATA TO A JSON FILE
    # #######################################  
    
    # if filename is not None:

    #     with open(filename, 'w') as f:
    #         json.dump(x, f,indent = 4)

    # return x












