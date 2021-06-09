# Redo cost function
# Redo cost differential function
# Redo backward solve
# Construct initial guess

from  numpy import *
from sysParamFns import *
from scipy.optimize import root,minimize_scalar
from scipy.integrate import solve_ivp
from sys import exit



def forwardSolveV2(xI,jt,nt,uTmp,dt,a,b,c,d,e,f,g,params,at,rt):
    
    # This routine does a forward solution, starting at time index jt up to nt
    # The u vector here contains the control information for this time interval
    nt=int(nt)
    jt = int(jt)
    xTmp=zeros((nt-jt+1,len(xI)))# functions of time
    xTmp[0,:] = xI
    xTmpH=zeros((nt-jt,len(xI)))#half-step solution needed for Simpson's rule integration    
 
    # Forward calculation of x, starting from time jt
    for j in range(nt-jt):
        tInt = (j+jt + arange(3)/2)*dt
        paramsTmp = params + [uTmp[j,:]]
        sol = solve_ivp(Afn, [tInt[0],tInt[-1]], xTmp[j,:], t_eval=tInt, args = (paramsTmp,),\
                              first_step=dt/2, atol=at, rtol = rt)
        tmp = shape(sol.y)# Check to make sure solution is complete
        if tmp[1]<3:
            # Improve tolerance to make sure a solution is found
            sol = solve_ivp(Afn, [tInt[0],tInt[-1]], xTmp[j,:], t_eval=tInt, args = (paramsTmp,),\
                            first_step=dt/2, atol=1, rtol = .001)
        if tmp[1]<3:
            cost = 1E100 # Solution has failed--bail out        
            return(xTmp,xTmpH,cost)
        else:    
            xInt = maximum(sol.y,0)# Ensure that solution is nonnegative
        
            xTmpH[j,:]=xInt[:,1]
            xTmp[j+1,:]=xInt[:,-1]
    cost,allCostsSummed = costFnV2(xI,xTmp,xTmpH,uTmp,dt,a,b,c,d,e,f,g,params)
    return xTmp,xTmpH,cost,allCostsSummed

    
def costFnV2(xI,xTmp,xTmpH,uTmp,dt,a,b,c,d,e,f,g,params):

    wY,wA,wPY,wPA,beta,sig,tau,rhoA,rhoY,gamA,gamY,gamH,Pi,eta,nu,mu,theta,rr,Phi = params 

    # Use Simpson's rule, which is used because of discontinuity in u.
    # Multiply by dt at the end
    # Notice that u[j] is the control on the interval [j*dt, (j+1)*dt]
    nt = size(xTmpH,axis=0)
    # Initialize arrays
    testCost = zeros((nt,2))
    distCost = zeros((nt,2))
    sickCost = zeros((nt,2))
    hospCost = zeros((nt,2))
    
    # Fixed testing cost subgroup 0
    testCost[:,0] = a[0,0]*(uTmp[:,0]>0) 

    # Linear and quadratic testing cost subgroup 0    
    sumTmp = sum(xTmp[:,:5],axis=1)
    uN_L = uTmp[:,0]*sumTmp[:-1]
    uN_H = uTmp[:,0]*sum(xTmpH[:,:5],axis=1)
    uN_R = uTmp[:,0]*sumTmp[1:]
    testCost[:,0] += a[0,1]*(uN_L + 4*uN_H + uN_R)/6 #Linear
    testCost[:,0] += a[0,2]*uTmp[:,0]*(uN_L + 4*uN_H  + uN_R)/6 #Quadratic
    
    # Fixed testing cost subgroup 1
    testCost[:,1] = a[1,0]*(uTmp[:,1]>0) 
    
    # Linear and quadratic testing cost subgroup 1    
    sumTmp = sum(xTmp[:,9:14],axis=1)
    uN_L = uTmp[:,1]*sumTmp[:-1]
    uN_H = uTmp[:,1]*sum(xTmpH[:,9:14],axis=1)
    uN_R = uTmp[:,1]*sumTmp[1:]
    testCost[:,1] += a[1,1]*(uN_L + 4*uN_H + uN_R)/6
    testCost[:,1] += a[1,2]*uTmp[:,1]*(uN_L + 4*uN_H  + uN_R)/6


    # Linear and quadratic distancing cost subgroup 0    
    sumTmp = sum(xTmp[:,:8],axis=1)
    uN_L = uTmp[:,2]*sumTmp[:-1]
    uN_H = uTmp[:,2]*sum(xTmpH[:,:8],axis=1)
    uN_R = uTmp[:,2]*sumTmp[1:]
    distCost[:,0] = b[0,1]*(uN_L + 4*uN_H + uN_R)/6 #Linear
    distCost[:,0] += b[0,2]*uTmp[:,2]*(uN_L + 4*uN_H  + uN_R)/6 #Quadratic
    
    # Linear and quadratic testing cost subgroup 1    
    sumTmp = sum(xTmp[:,9:-1],axis=1)
    uN_L = uTmp[:,3]*sumTmp[:-1]
    uN_H = uTmp[:,3]*sum(xTmpH[:,9:-1],axis=1)
    uN_R = uTmp[:,3]*sumTmp[1:]
    distCost[:,1] = b[1,1]*(uN_L + 4*uN_H + uN_R)/6
    distCost[:,1] += b[1,2]*uTmp[:,3]*(uN_L + 4*uN_H  + uN_R)/6

    
    sickCost[:,0] = c[0]*(xTmp[:-1,5]+xTmp[1:,5]+4*xTmpH[:,5])/6
    sickCost[:,1] = c[1]*(xTmp[:-1,14]+xTmp[1:,14]+4*xTmpH[:,14])/6

    hospCost[:,0] = d[0]*(xTmp[:-1,6]+xTmp[1:,6]+4*xTmpH[:,6])/6
    hospCost[:,1] = d[1]*(xTmp[:-1,15]+xTmp[1:,15]+4*xTmpH[:,15])/6

    # Summed costs per time interval from all sources    
    totCost = sum(testCost + distCost + sickCost + hospCost, axis=1)*dt 
    # Append the costs from total deaths (e) & remaining infected (f)

    # totCost = append(totCost,   e[0]*xTmp[-1,8] + e[1]*xTmp[-1,17] \
    #                           + f[0]*xTmp[-1,4] + f[1]*xTmp[-1,13] \
    #                           + g[0]*xTmp[-1,0] + g[1]*xTmp[-1,9])
    
    # totCost = append(totCost,   e[0]*xTmp[-1,8] + e[1]*xTmp[-1,17] \
    #                           + f[0]*xTmp[-1,4] + f[1]*xTmp[-1,13] \
    #                           + g[0]*( sum(xTmp[-1,0:7], axis = 0) ) + g[1]*( sum(xTmp[-1,9:16],axis = 0) ) )

    # Compartments: xI
    S =  [xTmp[-1,0], xTmp[-1,9] ]   # [0] Susceptible(S), 
    E =  [xTmp[-1,1], xTmp[-1,10] ]  # [1] Exposed (E), 
    PY = [xTmp[-1,2], xTmp[-1,11] ]  # [2] Pre-symptomatic infectious (PY ), 
    PA = [xTmp[-1,3], xTmp[-1,12] ]  # [3] Pre-asymptomatic infectious (PA), 
    IY = [xTmp[-1,4], xTmp[-1,13] ]  # [4] Symptomatic infectious (IY ), 
    IA = [xTmp[-1,5], xTmp[-1,14] ]  # [5] Asymptomatic infectious (IA), 
    IH = [xTmp[-1,6], xTmp[-1,15] ]  # [6] Symptomatic infectious that are hospitalized (IH), 
    R =  [xTmp[-1,7], xTmp[-1,16] ]  # [7] Recovered (R), 
    D =  [xTmp[-1,8], xTmp[-1,17] ]  # [8] Deceased (D) 

    # MAKE A V_prime ( v = nu )
    V_prime = (nu * mu) / ( (1-nu)*gamH + nu*mu )  

    # MAKE A Pi_prime
    Pi_prime = (Pi*eta) / ( (Pi * eta) + (1 - Pi) * gamY )

    # Costs at final time from total deaths and herd immunity situation
    # Low risk final costs
    C0herd = g[0]*( S[0] + f[0]*( (1-tau)*E[0] + PA[0] + IA[0] ) )
    C0death = e[0]*( D[0] + V_prime[0]*( IH[0] + Pi_prime[0]*( IY[0] + PY[0] + tau*E[0] ) ) )
    
    # High risk final costs
    C1herd = g[1]*( S[1] + f[1]*( (1-tau)*E[1] + PA[1] + IA[1] ) )
    C1death =  e[1]*( D[1] + V_prime[1]*( IH[1] + Pi_prime[1]*( IY[1] + PY[1] + tau*E[1] ) ) )

    totCost = append(totCost, C0herd + C0death + C1herd + C1death)


    # [1] Testing  (Low Risk)
    # [2] Testing (High Risk)
    # [3] Social Distance (Low Risk)
    # [4] Social Distance (High Risk)
    # [5] Opportunity of Sickness
    # [6] Hospital    
    # [7] Deaths (at tf)
    # [8] Remaining Infected (at tf)
    # [9] Remaining Susceptible (at tf)
    a1 = sum(testCost[:,0])*dt 
    a2 = sum(testCost[:,1])*dt 
    a3 = sum(distCost[:,0])*dt 
    a4 = sum(distCost[:,1])*dt 
    a5 = sum(sickCost)*dt 
    a6 = sum(hospCost)*dt 
    
    allCostsSummed = [a1,a2,a3,a4,a5,a6,C0herd+C1herd,C0death+C1death]
                        
        
    return totCost,allCostsSummed

    

def forwardSolve(xI,jt,nt,uTmp,dt,a,b,c,d,e,f,params,at,rt):
    
    # This routine does a forward solution, starting at time index jt up to nt
    # The u vector here contains the control information for this time interval
    nt=int(nt)
    jt = int(jt)
    xTmp=zeros((nt-jt+1,len(xI)))# functions of time
    xTmp[0,:] = xI
    xTmpH=zeros((nt-jt,len(xI)))#half-step solution needed for Simpson's rule integration    
 
    # Forward calculation of x, starting from time jt
    for j in range(nt-jt):
        tInt = (j+jt + arange(3)/2)*dt
        paramsTmp = params + [uTmp[j,:]]
        sol = solve_ivp(Afn, [tInt[0],tInt[-1]], xTmp[j,:], t_eval=tInt, args = (paramsTmp,),\
                              first_step=dt/2, atol=at, rtol = rt)
        tmp = shape(sol.y)# Check to make sure solution is complete
        if tmp[1]<3:
            # Improve tolerance to make sure a solution is found
            sol = solve_ivp(Afn, [tInt[0],tInt[-1]], xTmp[j,:], t_eval=tInt, args = (paramsTmp,),\
                            first_step=dt/2, atol=1, rtol = .001)
        if tmp[1]<3:
            cost = 1E100 # Solution has failed--bail out        
            return(xTmp,xTmpH,cost)
        else:    
            xInt = maximum(sol.y,0)# Ensure that solution is nonnegative
        
            xTmpH[j,:]=xInt[:,1]
            xTmp[j+1,:]=xInt[:,-1]
    cost,allCostsSummed = costFn(xI,xTmp,xTmpH,uTmp,dt,a,b,c,d,e,f)
    return xTmp,xTmpH,cost,allCostsSummed

def costFn(xI,xTmp,xTmpH,uTmp,dt,a,b,c,d,e,f):
    # Use Simpson's rule, which is used because of discontinuity in u.
    # Multiply by dt at the end
    # Notice that u[j] is the control on the interval [j*dt, (j+1)*dt]
    nt = size(xTmpH,axis=0)
    # Initialize arrays
    testCost = zeros((nt,2))
    distCost = zeros((nt,2))
    sickCost = zeros((nt,2))
    hospCost = zeros((nt,2))
    
    # Fixed testing cost subgroup 0
    testCost[:,0] = a[0,0]*(uTmp[:,0]>0) 

    # Linear and quadratic testing cost subgroup 0    
    sumTmp = sum(xTmp[:,:5],axis=1)
    uN_L = uTmp[:,0]*sumTmp[:-1]
    uN_H = uTmp[:,0]*sum(xTmpH[:,:5],axis=1)
    uN_R = uTmp[:,0]*sumTmp[1:]
    testCost[:,0] += a[0,1]*(uN_L + 4*uN_H + uN_R)/6 #Linear
    testCost[:,0] += a[0,2]*uTmp[:,0]*(uN_L + 4*uN_H  + uN_R)/6 #Quadratic
    
    # Fixed testing cost subgroup 1
    testCost[:,1] = a[1,0]*(uTmp[:,1]>0) 
    
    # Linear and quadratic testing cost subgroup 1    
    sumTmp = sum(xTmp[:,9:14],axis=1)
    uN_L = uTmp[:,1]*sumTmp[:-1]
    uN_H = uTmp[:,1]*sum(xTmpH[:,9:14],axis=1)
    uN_R = uTmp[:,1]*sumTmp[1:]
    testCost[:,1] += a[1,1]*(uN_L + 4*uN_H + uN_R)/6
    testCost[:,1] += a[1,2]*uTmp[:,1]*(uN_L + 4*uN_H  + uN_R)/6


    # Linear and quadratic distancing cost subgroup 0    
    sumTmp = sum(xTmp[:,:8],axis=1)
    uN_L = uTmp[:,2]*sumTmp[:-1]
    uN_H = uTmp[:,2]*sum(xTmpH[:,:8],axis=1)
    uN_R = uTmp[:,2]*sumTmp[1:]
    distCost[:,0] = b[0,1]*(uN_L + 4*uN_H + uN_R)/6 #Linear
    distCost[:,0] += b[0,2]*uTmp[:,2]*(uN_L + 4*uN_H  + uN_R)/6 #Quadratic
    
    # Linear and quadratic testing cost subgroup 1    
    sumTmp = sum(xTmp[:,9:-1],axis=1)
    uN_L = uTmp[:,3]*sumTmp[:-1]
    uN_H = uTmp[:,3]*sum(xTmpH[:,9:-1],axis=1)
    uN_R = uTmp[:,3]*sumTmp[1:]
    distCost[:,1] = b[1,1]*(uN_L + 4*uN_H + uN_R)/6
    distCost[:,1] += b[1,2]*uTmp[:,3]*(uN_L + 4*uN_H  + uN_R)/6

    
    sickCost[:,0] = c[0]*(xTmp[:-1,5]+xTmp[1:,5]+4*xTmpH[:,5])/6
    sickCost[:,1] = c[1]*(xTmp[:-1,14]+xTmp[1:,14]+4*xTmpH[:,14])/6

    hospCost[:,0] = d[0]*(xTmp[:-1,6]+xTmp[1:,6]+4*xTmpH[:,6])/6
    hospCost[:,1] = d[1]*(xTmp[:-1,15]+xTmp[1:,15]+4*xTmpH[:,15])/6

    # Summed costs per time interval from all sources    
    totCost = sum(testCost + distCost + sickCost + hospCost, axis=1)*dt 
    # Append the costs from total deaths (e) & remaining infected (f)
    totCost = append(totCost,   e[0]*(xTmp[-1,8]-xI[8]) + e[1]*(xTmp[-1,17]-xI[17]) \
                              + f[0]*(xTmp[-1,4]-xI[4]) + f[1]*(xTmp[-1,13]-xI[13]) )
    
    # [1] Testing  (Low Risk)
    # [2] Testing (High Risk)
    # [3] Social Distance (Low Risk)
    # [4] Social Distance (High Risk)
    # [5] Opportunity of Sickness
    # [6] Hospital
    # [7] Deaths (during interval)
    # [8] Remaining Infected (additonal incurred)
    allCostsSummed = [sum(testCost[:,0])*dt,sum(testCost[:,1])*dt,\
                      sum(distCost[:,0])*dt,sum(distCost[:,1])*dt,\
                      sum(sickCost)*dt, sum(hospCost)*dt,\
                      e[0]*(xTmp[-1,8]-xI[8])+e[1]*(xTmp[-1,17]-xI[17]),\
                      f[0]*(xTmp[-1,4]-xI[4])+ f[1]*(xTmp[-1,13]-xI[13]) ]
                        
        
    return totCost,allCostsSummed
    
