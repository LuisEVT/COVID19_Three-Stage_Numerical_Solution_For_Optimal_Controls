# Use one interval for enum
# Modify cost function & let CT verify
# Save data to produce MC improvement graphs


# To document: update of probability for pairs of intervals
# Update of probability for treatments for chosen intervals
# Temperature adjustment
# 

### LOCAL IMPORTS
from controlfnsenumv5 import forwardSolveV2
from initial import initial 

### OTHER IMPORTS
import numpy as np
from skimage import measure # USED FOR CONTOUR LINES
from scipy.sparse.linalg import eigs as speigs

import time
import matplotlib.pyplot as plt


# This algorithm will randomly choose a pair of intervals, 
# and randomly perturb the controls on that pair of intervals.
# The algorithm is designed so as to be flexible, as well as 
# relatively few parameters to optimize.
def monteCarlo(a,b,c,d,e,f,g,params,tf,xI,u,uMax,nInt = 20,tempParam = {'flatIter':10,'NoImproveMax':3} ,at= 10,rt= 0.001,initStep = 0.1):
    '''
    Parameters:     
        nInt: Number of intervals on which MC changes are made
        u: Current trial control values    
        uMax: Maximum control values
        xI: Initial values
        tf: Final time
        a: Test Cost
        b: Distance Cost
        c: Sickness Cost
        d: Hospital Cost
        e: Death Cost
        f: Residual Infection Cost
        params: system parameters for the solution
        at: Absolute Tolerance in solution
        rt: Relative Tolerance in solution
        initStep: Initial change in control, if t
        nIter:  Number of Monte Carlo Iterations
    '''
    # OLD NUMBER OF STEPS INTERVALS & NUMBER OF CONTROLS
    nStep , nCtrl = u.shape 


    # NUMBER OF STEPS PER INTERVAL
    # NOTE NUMBER OF INTERVALS SHOULD DIVIDE EVENLY
    intStep = tf//nInt
    if intStep != tf/nInt:
        print("WARNING-MC INTERVALS DO NOT DIVIDE THE TIME SCALE EVENLY")

    # INTEGRATION STEP
    dt = tf/nStep # day step

    startStep = 0 # starting interval for recomputing cost

    #################

    # BEST CONTROL MEASURES PER INTERVAL 
    # ( FOR EACH INTERVAL SET, AVERAGE THE CONTROL MEASURES )
    uBestInt = np.array([ np.mean( u[ ii:ii+intStep,: ] , axis = 0)   for ii in range(0,tf,intStep)])
    
    # BEST CONTROL MEASURES PER INT EXPANDED TO INT STEP
    uBestStep = np.repeat(uBestInt,intStep,axis=0)
    
    # GET MY INITIAL SOLUTION TO WORK WITH
    x, xH, costPerStep, _ = forwardSolveV2(xI,0,tf,uBestStep,dt,a,b,c,d,e,f,g,params,at,rt)

    # INITIAL BEST COST
    costBest = sum(costPerStep)

    ###########
    # INITIALIZE PROBABILITY VECTORS 
    ###########

    # INTERVAL PAIR SELECTION
    pPair0 = np.ones( (nInt,nInt) ) + np.eye(nInt) 
    pPair = 1.*pPair0
    
    # CONTROL CHANGE PROBABILITY BY INTERVAL    
    pChange0 = 0.5 * np.ones( (nInt,nCtrl) ) 
    pChange = 1.*pChange0
    
    # MAXIMUM SIZE OF RANDOM CONTROL CHANGE
    randStep = np.zeros( (nInt,nCtrl) ) 

    pPairMin = 0.5
    deltaP = 0.5
    deltaPair = 0.5
    deltaStep = 0.5

    # Temperature parameters
    # WHEN CURFLATITER = FLATITER, THEN INCREASE TEMPERATURE
    #flatIter = 10 #25 #30
    flatIter = tempParam['flatIter']
    # Number of flat intervals with no improvement triggers end of MC
    #noImproveMax = 3 #3 #5 # 
    noImproveMax = tempParam['noImproveMax']
    # Memory factor used in computing the mean increase
    meanMemoryFactor = 0.95
    # Temperature scale factor (multiplied by estimated mean increase)
    tempMultiplierReset = 4
    # Temperature incremental mult. factor
    tempReduceFactor = 0.95

    # INITIAL VALUES
    curFlatIter = 0 # CURRENT NUMBER OF PREVIOUS ITERATIONS WITH NO DECREASE
    tempMultiplier = 0.001 #  CLOSE TO BUT NOT = 0, TO AVOID INVALID VALUES
    meanIncrEst = 0.001 #  CLOSE TO BUT NOT = 0, TO AVOID INVALID VALUES
    
    uInt = 1.0*uBestInt  # CONTROL MEASURES PER INTERVAL

    costBestArr = [costBest]
    iterCost = [costBest]

    kk = 0
    noImproveCount = 0    
    
    # Set best previous cost
    costBestPrev = np.inf
    
    # count number of temperature increases that don't bring any improvement
    # to overall best.  When reach limit, then stop
    while (noImproveCount < noImproveMax):

        ####################################
        # CHANGE CONTROLS ON INTERVAL IX,JX
        ####################################

        # CONTROL CHANGES WILL BE RECORDED HERE
        randStepVec = np.zeros( (nCtrl,2) ) # control changes will be recorded here
        chVec = 0.*randStepVec
        #t1 = time.perf_counter()
        
        # MAKE SURE THAT AT LEAST ONE OF THE CONTROLS HAS CHANGED IN EACH INTERVAL
        while ( (np.amax(chVec[:,0])) + (np.amax(chVec[:,1])) == 0):
        
            # SELECT 2 RANDOM INTERVALS FOR CHANGE
            pPairVec = pPair.flatten()
            pPairVec = pPairVec / np.sum(pPair)
            randIx =  np.random.choice( range(len(pPairVec))  , size = 1 , p = pPairVec )[0]
    
            # TWO INTERVALS SELECTED FOR CHANGE
            ix = [randIx % nInt, randIx // nInt] 
    
            uNew = 1.0 * uInt
    
            for ii in range(2):
                
                # DO THE SAME THING FOR J
                if (ii==0) or not(ix[0] == ix[1]):
                    
                    # CHOOSE WHICH CONTROLS WILL CHANGE
                    # Normalize sum of probabilities to 1
                    # (This is to make sure probabilities don't get too small)
                    pChange[ix[ii],:] /= sum(pChange[ix[ii],:])
                    chVec[:,ii] = np.random.binomial(n = 1, p = pChange[ix[ii],:]) 
    
                    # SIGN OF CHANGE
                    # 0.5 is because comparing to previous step.
                    stepSignVec = np.random.choice([-0.5,1], size = 4, p = [1/3,2/3]) 
            
                    # SIZE OF STEP (TO BE SCALED)
                    unifVec = np.random.uniform(size = 4)
    
                    # ACTUAL STEP SIZE
                    # IF YOU START AT ZERO CONTROL MEASURES, START OFF WITH SOMETHING
                    randStepVec[:,ii] = chVec[:,ii] * randStep[ix[ii],:] * unifVec * stepSignVec \
                        + chVec[:,ii]*( randStep[ix[ii],:] == 0 ) * unifVec * np.sign(stepSignVec) * initStep
            
                    uNew[ix[ii],:] = uNew[ix[ii],:] + randStepVec[:,ii]
                    uNew[ix[ii],:] = np.minimum(uNew[ix[ii],:] , uMax)
                    uNew[ix[ii],:] = np.maximum(uNew[ix[ii],:] , np.zeros(nCtrl))


        
                    # print(uNew[ix[ii],:])
                    # t = 1
        # print('\n')
        #t2 = time.perf_counter()

        ####################################
        # COMPUTE COST USING FORWARD 
        ####################################       
        startInt = min(ix) 
        # Recompute solution from here.
        startStep = startInt * intStep

        uTmp = np.copy(np.repeat( uNew, intStep, axis = 0 ))

        xTmp,xTmpH,costPerStepTmp,_ = forwardSolveV2( x[startStep,:], \
                                                  startStep, nStep, \
                                                  uTmp[startStep:,:], \
                                                  dt,a,b,c,d,e,f,g,params,at,rt)

        #t3 = time.perf_counter()

        # EVALUATE CHANGE IN COST
        costDiff = sum(costPerStepTmp)  - sum(costPerStep[startStep:])

        # SAVE THE CURRENT COST IN THIS ITERATION
        iterCost.append( sum(costPerStepTmp)  + sum(costPerStep[:startStep]) )

        # print('\nIter:{} CostDiff:{:.4E} tempM:{:.4E} flat:{}/{:.1f}'.format(kk,costDiff,tempMultiplier,curFlatIter,flatIter))

        # DECIDE IF NEW STRATEGY IS ACCEPTED
        if costDiff < 0: # ACCEPT CHANGE, COST HAS IMPROVED

            # RESET NO-IMPROVEMENT COUNT            
            curFlatIter = 0

            ## UPDATE SOLUTION, CONTROL & COST
            
            # UPDATE CONTROL VALUES ON INTERVALS
            uInt = 1.0*uNew
            # Update solution values
            x[startStep:,:] = 1.0*xTmp
            xH[startStep:,:] = 1.0*xTmpH
            # UPDATE COST
            costPerStep[startStep:] = 1.0*costPerStepTmp

            # UPDATE PROBABILITY OF CHANGING CONTROL 'pChange'
            # INCREASED PROBABILITY OF 'GOOD' CHOICE
            pChange[ix[0],:] += deltaP * chVec[:,0]*(1-pChange[ix[0],:]) 
            
            # REPEAT SUCCESSFUL STEP AND INCREASE STEP
            randStep[ix[0],:] = randStepVec[:,0] / unifVec * (1 + deltaStep) 

            if not(ix[0] == ix[1]): # IF INTERVALS ARE DIFFERENT, MAKE SAME CHANGE FOR 2ND INTERVAL
                pChange[ix[1],:] += deltaP * chVec[:,1]*(1-pChange[ix[1],:])
                randStep[ix[1],:] = randStepVec[:,1] / unifVec * (1 + deltaStep)


            # UPDATE 'pPair'
            pPair[ix[0],ix[1]] += deltaPair

            #### @@@ ADD CONDITION WHEN COST IS BETTER THAN BEST PREVIOUS ACHIEVED
            #### UPDATE UBEST AND BEST COST
            if sum(costPerStep) < costBest:
                uBestInt = 1.0*uInt
                costBest = sum(costPerStep)
            
        else: # NEW CONTROL DOES NOT LOWER COST
      
            # REDUCE PROBABILITY OF CHANGING CONTROL 'pChange'
            pChange[ix[0],:] -= deltaP * chVec[:,0] * pChange[ix[0],:] # DECREASED PROBABILITY OF 'GOOD' CHOICE
            
            if not(ix[0] == ix[1]):
                pChange[ix[1],:] -= deltaP * chVec[:,1] * pChange[ix[1],:] 


            # REDUCE PROBABILITY OF CHOOSING 'pPair'
            pPair[ix[0],ix[1]] -= deltaPair
            pPair[ix[0],ix[1]] = max( pPair[ ix[0],ix[1] ] , pPairMin)

            # INCREASE COUNT OF NON-DECREASING ITERATIONS
            curFlatIter += 1
            # UPDATE ESTIMATE OF MEAN INCREASE
            meanIncrEst = meanMemoryFactor * meanIncrEst + (1 - meanMemoryFactor) * costDiff

            # CHECK IF IMPROVEMENTS IN COST HAS BEEN STATIC FOR A WHILE
            if curFlatIter >= flatIter :
                # RESET PROBABILITIES
                pPair = 1.*pPair0
                pChange = 1.*pChange0
                # INCREASE TEMPERATURE TO ENCOURAGE EXPLORATION 
                tempMultiplier = tempMultiplierReset
                # CHANGE NON-INCREASE CONSECUTIVE ITERATIONS TOLERANCE
                # flatIter *= flatIterChangeFactor
                # RESET FLAT ITERATION COUNT
                curFlatIter = 0
                
                #See if this temp cycle brought improvement
                #If improvement, then reset noImproveCount
                if costBest < costBestPrev:
                    noImproveCount = 0 #max(noImproveCount - 1,0)
                    costBestPrev = costBest
                #If no improvement, then increment noImproveCount
                #MC will stop when noImproveCount reaches threshold
                else:
                    noImproveCount +=1


            # ACCEPT CHANGE ANYWAY, DEPENDING ON TEMPERATURE
            if (np.random.uniform() < np.exp( -costDiff / (meanIncrEst*tempMultiplier))):
                
                # Update solution, control & cost
                # Update control values on intervals
                uInt = 1.0*uNew
                # Update solution values
                x[startStep:,:] = 1.0*xTmp
                xH[startStep:,:] = 1.0*xTmpH
                # Update cost
                costPerStep[startStep:] = 1.0*costPerStepTmp

        #t4 = time.perf_counter()
        #print('CCT:{:.4E}  FST:{:.4E} IT:{:.4E}'.format(t2-t1,t3-t2,t4-t1))

        # PRINT FOR DEBUG
        #print('Iter:{} Cost:{:.8f} flatCount:{} NoImproveCount:{}'.format(kk,sum(costPerStep),curFlatIter,noImproveCount))
        
        # REDUCE THE TEMPERATURE 
        tempMultiplier *= tempReduceFactor
        kk += 1
        costBestArr.append(costBest)
        

    return np.repeat(uBestInt,intStep,axis=0) , costBest , costBestArr , iterCost


def enumOpt(a,b,c,d,e,f,g,params,tf,nt,xI,nsPerInt,ctrlVals,ctrlEq,uMax,at= 10,rt= 0.001,**kwargs):
    '''
    Returns:
        1. uBest: Controls values on integration steps ( number of integration steps, number of controls )
        2. costBest: optimal minimum cost
    Parameters:
        a: Test Cost
        b: Distance Cost
        c: Sickness Cost
        d: Hospital Cost
        e: Death Cost
        f: Residual Infection Cost
        params: List of System Parameters
        tf: Final time
        nt: Number of time intervals
        xI: Initial population with given characteristics ( low-risk & high-risk populations)
        nsPerInt: Number of integration steps per interval
        ctrlVals:  Possible control values in enumerative checking
        ctrlEq: Which controls are equal to others
        uMax: Maximum values for control measures 
        at: Absolute Tolerance in solution
        rt: Relative Tolerance in solution

    Defined in subroutine:
        uBest: Initial control values on integration steps
        costBest: Current optimal total cost
    '''

 
    # Initialize
    # Parameters for enumerative
    uBest = np.zeros((nt,4))  # initial control values on integration steps
    costBest = 1E100        # Set large initial value  

    
    # NUMBER OF TIME STEP , NUMBER OF CONTROLS
    nStep,nCtr = np.shape(uBest) 

    # TOTAL NUMBER OF INTERVAL IN ENUMERATIVE
    nEnumInt = len(nsPerInt)

    # TOTAL NUMBER OF TIME STEPS (SHOULD MATCH 'nt')
    tmp = np.cumsum(nsPerInt)
    if tmp[-1] != nStep:
        print("Error! nStep value doesn't agree with nsPerInt")   

    # GIVE START AND FINISH STEP INDEX FOR EACH INTERVAL
    isx = np.zeros(nEnumInt+1 , dtype = int)
    isx[1:] = tmp
    
    dt = tf/nStep

    #########
    # Set up vector to store control at integration step resolution
    ######### 

    # A COPY OF CONTROLS FOR EACH TIME STEP ( full-length control )
    uTmp = 1.0*uBest 

    # GIVES THE NUMBER OF TRIAL CONTROLS OF EACH TYPE ON EACH INTERVAL
    nCtrl = np.zeros((nEnumInt,nCtr), dtype = int)

    # Total number of different controls tried during enumerative alg.
    nCtrlTot = 1

    for ji in range(nEnumInt):
        # SET VARIABLE CONTROLS ON INTERVAL
        for ju in range (nCtr):
            # NUMBER OF TRIAL CONTROLS
            nCtrl[ji,ju] = len(ctrlVals[ji][ju])
            # ASSIGN VALUE            
            if nCtrl[ji,ju] > 0:
                nCtrlTot *= nCtrl[ji,ju]
                uTmp[isx[ji]:isx[ji+1],ju] = ctrlVals[ji][ju][0]
        for ju in range (nCtr):
            # NUMBER OF TRIAL CONTROLS
            nCtrl[ji,ju] = len(ctrlVals[ji][ju])
            # ASSIGN SPECIFIC VALUES
            if nCtrl[ji,ju] == 0:
                uTmp[isx[ji]:isx[ji+1],ju] = ctrlVals[ji][ctrlEq[ji][ju]][0]
                

    # COSTS ARE EVALUATED PER INTERVAL    
    # MOST RECENT COST PER INTERVAL
    costPerStep = np.zeros(nStep+1)
       
    # THIS IS THE SOLUTION
    x = np.zeros((nt+1,len(xI)))
    xH = np.zeros((nt,len(xI)))
    x[0,:] = xI

        
    # TRY ALL CONTROL COMBINATIONS, AND KEEP THE BEST
    for jc in range(nCtrlTot):
        # print(jc)
        # construct initial values for each control on each interval
        # start from the back
        jcTmp = jc
        jcTmp1 = (jc-1)%nCtrlTot # shows which interval/control comb. changes

        # Update intervals backward so that changes are towards the end
        for ji in range(nEnumInt-1,-1,-1):
            # new starting step
            startStep = int(sum(nsPerInt[:ji]))
            # Update the controls: look at all controls in interval
            for ju in range(nCtr):
                # Check and see if control is variable
                if nCtrl[ji,ju] > 0:
                    # Set control equal to next value
                    uIntTmp = ctrlVals[ji][ju][jcTmp % nCtrl[ji,ju] ]
                    # Reset value of control on interval
                    uTmp[startStep:startStep+nsPerInt[ji],ju] = uIntTmp
                    # Update change index
                    jcTmp = jcTmp //nCtrl[ji,ju]
                    jcTmp1 = jcTmp1 // nCtrl[ji,ju]
                    # Check if finished
                    if jcTmp == jcTmp1:
                        break
            for ju in range(nCtr):
                # Check and see if control is variable
                if nCtrl[ji,ju] == 0:
                    uTmp[startStep:startStep+nsPerInt[ji],ju] = uTmp[startStep,ctrlEq[ji][ju]]
            # If finished, then terminate
            if jcTmp == jcTmp1:
                break


        xTmp,xTmpH,costPerStepTmp,_ = forwardSolveV2( x[startStep,:], \
                                                      startStep, nStep, \
                                                      uTmp[startStep:,:], \
                                                      dt,a,b,c,d,e,f,g,params,at,rt)
            
        # UPDATE COST PER STEP    
        costPerStep[startStep:] = 1.0*costPerStepTmp

        # UPDATE SOLUTION (x and xH)
        x[startStep:,:] = 1.0*xTmp
        xH[startStep:,:] = 1.0*xTmpH        
        

        # EVALUATE CHANGE IN COST
        # NOTE cost function does not include deaths and remaining
        # infected at start, so they must be added in.
        costDiff = sum(costPerStep)  - costBest
         
        # Decide if new strategy is accepted
        if costDiff <= 0: # Accept change, cost has improved

            # Update control values on intervals
            uBest = 1.0*uTmp

            # LV: VERIFY THIS IS CORRECT
            costBest = sum(costPerStep)

    # Return best control per integration step   
    return uBest,costBest


def localOpt(a,b,c,d,e,f,g,params,tf,nt,uBest,uMax,xI,du = 0.01,utol = 0.01 ,at= 10,rt= 0.001, nIterMax = 200,**kwargs):

    # OLD NUMBER OF STEPS INTERVALS & NUMBER OF CONTROLS
    nStep , nCtrl = uBest.shape 

    # INTEGRATION STEP
    dt = tf/nStep # day step

    # COMPUTE INITIAL SOLUTION
    xBest,xBestH,costPerStepBest,_ = \
        forwardSolveV2(xI,0,tf,uBest,dt,a,b,c,d,e,f,g,params,at,rt)

    allPts = False    
    for ji in range(nIterMax): # Iterate until convergence
        #print(ji,':',nIterMax)
        #print('iter:{} out of {}'.format(ji+1,nIterMax))
        for iter in range(2):  # Forward and backward
            if iter==0:
                jtList = np.arange(nt)
            else:
                jtList = np.arange(nt,0,-1)-1## go backwards
            
            # Set benchmark cost-- if no improvement, then quit
            costBench = sum(costPerStepBest)
                
            for jt in jtList:
                # Do the controls in random order
                oTmp = np.random.permutation(4)
                for jc in oTmp:
                    if allPts or (jt==0) or (jt==nt-1)\
                        or (jt > 0 and jt < nt-1 \
                            and (uBest[jt,jc] != uBest[jt-1,jc] \
                                  or uBest[jt,jc] != uBest[jt+1,jc])):
    
                        # Find locally optimal control
                        # Est. of 2nd derivative of objective function
                        # Four different controls
                        if jc==0: # Testing low risk
                            d2J = 2*dt*a[0,2]* \
                                  (sum(xBest[jt,:5]) + \
                                    4*sum(xBestH[jt,:5]) + \
                                    sum(xBest[jt+1,:5]))/6
                        elif jc==1: # Testing high risk
                            d2J = 2*dt*a[1,2]* \
                                    (sum(xBest[jt,9:14]) + \
                                    4*sum(xBestH[jt,9:14]) + \
                                    sum(xBest[jt+1,9:14]))/6
                        else: # distancing
                            d2J = 2*dt*b[jc-2,2]
    
                        # Est. of first derivative of objective fn.
                        # Perturb control                            
                        uTmp = np.copy(uBest[jt:,:])
                        uTmp[0,jc]+=du
                        # Compute objective with perturbed control
                        xTmp,xTmpH,costPerStepTmp,_ = \
                            forwardSolveV2(xBest[jt,:],jt,nt,uTmp,dt,a,b,c,d,e,f,g,params,at,rt)

                        costBestTmp = sum(costPerStepBest[jt:])
                        costDiff = sum(costPerStepTmp)- costBestTmp
                        # Correct for discontinuity at 0
                        if (jc<2 and uBest[jt,jc]==0):
                            costDiff += a[jc,0]*dt
                        # First deriv est. & correct for second derivative
                        d1J = costDiff/du - d2J*du/2

                        # New est. of locally  optimal control
                        if d2J == 0:
                            uNew = 0
                        else:
                            uNew = max( min( uBest[jt,jc]-d1J/d2J, np.max(uMax[jc]) ), 0 )


                        if abs(uBest[jt,jc] - uNew) > utol:
                            uTmp = np.copy(uBest[jt:,:])
                            uTmp[0,jc] = uNew
                            xTmp,xTmpH,costPerStepTmp,_ = \
                                forwardSolveV2(xBest[jt,:],jt,nt,uTmp,dt,a,b,c,d,e,f,g,params,at,rt)
                            
                            # If the new u actually reduces cost, then use it

                            if sum(costPerStepTmp) < costBestTmp:
                                uBest[jt,jc] = uNew
                                xBest[jt:,:]=xTmp
                                xBestH[jt:,:]=xTmpH
                                costPerStepBest[jt:]=costPerStepTmp

                        # (If necessary) compare with 0
                        if uBest[jt,jc] > 0 and jc<2:
                            uTmp = np.copy(uBest[jt:,:])
                            uTmp[0,jc]=0
                            xTmp,xTmpH,costPerStepTmp,_ = \
                                forwardSolveV2(xBest[jt,:],jt,nt,uTmp,dt,a,b,c,d,e,f,g,params,at,rt)

                            # If 0 control reduces cost, then use it

                            if sum(costPerStepTmp)<sum(costPerStepBest[jt:]):
                                uBest[jt,jc] = 0
                                xBest[jt:,:] = xTmp
                                xBestH[jt:,:] = xTmpH
                                costPerStepBest[jt:] = costPerStepTmp
                        
                        printout = [ji,jt,jc,uBest[jt,jc],sum(costPerStepBest)]
                        #print(*printout)
    
            # If no improvement, leave loop
            if costBench - sum(costPerStepBest) <0.0001*costBench:
                break
        # If no improvement, either check all points or leave
        if costBench - sum(costPerStepBest) <0.0001*costBench:
            if allPts:
                break
            else:
                allPts = True

    # End loop
    
    # evaluate once more to get costBest
    xTmp,xTmpH,costPerStepTmp,allCostsSummed = \
    forwardSolveV2(xI,0,tf,uBest,dt,a,b,c,d,e,f,g,params,at,rt)
    
    return uBest,xBest,costPerStepTmp,allCostsSummed
 
    
def optMinCost(a,b,c,d,e,f,g,params,tf,nt,xI,nsPerInt,cVals,ctrlVals,ctrlEq,uMax, enumOn = True, mcOn = True,**kwargs):
    '''
    Returns:

    Parameters:
        a: Test Cost
        b: Distance Cost
        c: Sickness Cost
        d: Hospital Cost
        e: Death Cost
        f: Residual Infection Cost
        params: List of System Parameters
        tf: Final time
        nt: Number of time intervals
        xI: Initial population with given characteristics ( low-risk & high-risk populations)
        nsPerInt: Number of integration steps per interval
        ctrlVals:  Possible control values in enumerative checking
        ctrlEq: Which controls are equal to others
        uMax: Maximum values for control measures 
        **kwargs:
            at: Absolute Tolerance in solution
            rt: Relative Tolerance in solution
            nInt: Number of intervals of Monte Carlo
            initStep: 
            nIterMax:
            du: Change in control used to estimate derivative
            utol: tolerance for u change

    Defined in subroutine:
        uBest: Initial control values on integration steps
        costBest: Current optimal total cost
    '''

    t1 = time.perf_counter()
    if enumOn:
        enumCtrl, enumCost = enumOpt(a,b,c,d,e,f,g,params,tf,nt,xI,nsPerInt,ctrlVals,ctrlEq,uMax,**kwargs)
    else:
        enumCtrl = np.zeros( (nt,4) )
        enumCost = 0
    t2 = time.perf_counter()

    #print('enumCost:{:.4E} enumTime:{:.4E} min'.format(enumCost, (t2-t1)/ 60) )

    if mcOn and enumOn:
        mcCtrl, mcCost , costBestArr , iterCost = monteCarlo(a,b,c,d,e,f,g,params,tf,xI,enumCtrl,uMax,**kwargs)
    elif mcOn:
        mcCtrl = np.zeros( (nt,4) )
        mcCtrl, mcCost , costBestArr , iterCost = monteCarlo(a,b,c,d,e,f,g,params,tf,xI,mcCtrl,uMax,**kwargs)
    else:
        mcCtrl = np.zeros( (nt,4) )
        mcCost = 0
        costBestArr = [0]
        iterCost = [0]

    #mcCtrl, mcCost , costBestArr , iterCost = [ 0*enumCtrl , 0 , [0] , [0] ]

    t3 = time.perf_counter()
    #print('mcCost:{:.4E} mcTime:{:.4E} min'.format(mcCost, (t3-t2) / 60 ) )

    if mcOn:
        optCtrl,xBest,costPerStep,allCostsSummed = localOpt(a,b,c,d,e,f,g,params,tf,nt,mcCtrl,uMax,xI,**kwargs)
        optCost = sum(costPerStep)
    elif enumOn and not(mcOn):
        optCtrl,xBest,costPerStep,allCostsSummed = localOpt(a,b,c,d,e,f,g,params,tf,nt,enumCtrl,uMax,xI,**kwargs)
        optCost = sum(costPerStep)    
    else:
        localCtrl = np.zeros( (tf,4) )
        optCtrl,xBest,costPerStep,allCostsSummed = localOpt(a,b,c,d,e,f,g,params,tf,nt,localCtrl,uMax,xI,**kwargs)
        optCost = sum(costPerStep)             

    t4 = time.perf_counter()

    runtimes = [ t2 - t1 , t3 - t2 , t4 - t3 ]


    #print('optCost:{:.4E} optTime:{:.4E} min'.format( sum(costPerStep), (t4-t3) / 60 ) )
    #print('total time:{:.4E}m'.format( (t4-t1) / 60 ) )

    # plt.figure(1)
    # plt.plot(range(len(costBestArr)),costBestArr , label = 'Min Cost')
    # plt.plot(range(len(costBestArr)),iterCost, label = 'Iter Cost')

    # #plt.plot(0,enumCost,'ro')

    # plt.xlabel('nIterations')
    # plt.ylabel('$USD Cost')
    # plt.legend()


    # plt.figure(2)
    # size , _ = mcCtrl.shape

    # plt.plot( range(size) , mcCtrl[:,0] , label = 'Low-Risk Testing')
    # plt.plot( range(size) , mcCtrl[:,1] , label = 'High-Risk Testing')
    # plt.plot( range(size) , mcCtrl[:,2] , label = 'Low-Risk Distancing')
    # plt.plot( range(size) , mcCtrl[:,3] , label = 'High-Risk Distancing')

    # plt.legend()
    # plt.title('MC Control Measures')
    # plt.xlabel('Days')
    # plt.ylabel('Control rate')

    
    # plt.figure(3)
    # size , _ = optCtrl.shape

    # plt.plot( range(size) , optCtrl[:,0] , label = 'Low-Risk Testing')
    # plt.plot( range(size) , optCtrl[:,1] , label = 'High-Risk Testing')
    # plt.plot( range(size) , optCtrl[:,2] , label = 'Low-Risk Distancing')
    # plt.plot( range(size) , optCtrl[:,3] , label = 'High-Risk Distancing')

    # plt.legend()
    # plt.title('Local Opt. Control Measures')
    # plt.xlabel('Days')
    # plt.ylabel('Control rate')

    # plt.show()
 
    return [enumCtrl, enumCost], [mcCtrl, mcCost], [optCtrl,optCost,xBest,allCostsSummed],[costBestArr , iterCost],runtimes



# CALL LOCAL OPTIMIZATION 
if __name__ == '__main__':
    
    # SINGLE RUN WITH GIVEN CONFIGURATIONS

    dictVar = initial()

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

    ###########
    ## NOTE: THE SIMULATIONS IS RUN IN Testbench.py
    ###########

    print('#####')
    print('Two Intervals')
    print('#####')
    # nsPerInt = [nt]
    nsPerInt = [nt//2, nt - (nt//2)]
    # nsPerInt = [nt//3, nt//3, nt - 2 * nt//3] 
    # nsPerInt = [nt//4, nt//4, nt//4, nt - 3 * nt//4]
    
    nIntervals = len(nsPerInt)

    # Enumerative: discrete values on discrete intervals
    ctrlVals = [ [ cVals[0,:], [], cVals[2,:], [] ]  for ii in range(nIntervals) ]
    # Which controls are equal to each other
    ctrlEq = [ [0,0,0,2]  for ii in range(nIntervals) ]      

    # Temperature Parameters
    tempParam = {'flatIter':10,
                 'noImproveMax':3}

    e = [10**6, 0.75*(10**6)]
    g = [10**4 , 10**4 ]
    optMinCost(a,b,c,d,e,f,g,params,tf,nt,xI,nsPerInt,cVals,ctrlVals,ctrlEq,uMax, 
    tempParam =tempParam, enumOn = True,mcOn = True)


    ## RUN 2
    print('#####')
    print('enum Three Intervals')
    print('#####')
    
    # nsPerInt = [nt]
    # nsPerInt = [nt//2, nt - (nt//2)]
    nsPerInt = [nt//3, nt//3, nt - 2 * nt//3] 
    # nsPerInt = [nt//4, nt//4, nt//4, nt - 3 * nt//4]
    
    nIntervals = len(nsPerInt)

    # Enumerative: discrete values on discrete intervals
    ctrlVals = [ [ cVals[0,:], [], cVals[2,:], [] ]  for ii in range(nIntervals) ]
    # Which controls are equal to each other
    ctrlEq = [ [0,0,0,2]  for ii in range(nIntervals) ]      

    # Temperature Parameters
    tempParam = {'flatIter':10,
                 'noImproveMax':3}

    e = [10**6, 0.75*(10**6)]
    g = [10**4 , 10**4 ]
    optMinCost(a,b,c,d,e,f,g,params,tf,nt,xI,nsPerInt,cVals,ctrlVals,ctrlEq,uMax, 
    tempParam =tempParam, enumOn = True,mcOn = True)
        