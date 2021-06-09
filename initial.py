import numpy as np


def initial():

    ##############
    # INITIAL CONDITIONS
    ##############

    N = [1340000,423000]
    NA = [1340000,423000]
    xI = np.zeros(18)
    xI[0] = N[0] # Susceptible low risk
    xI[9] = N[1] # Susceptible high risk
    xI[1] = 150
    xI[10] = 50
    xI[0] = N[0]-xI[1] # Susceptible low risk
    xI[9] = N[1]-xI[10] # Susceptible high risk
    
    #############
    # SYSTEM PARAMETERS
    #############

    # BASELINE TRANSMISSION RATE
    beta = 0.0640
    # CONTACT MATRIX
    Phi = np.array([[10.56,2.77],
                    [9.4,2.63]]) 
    # REVOVERY RATE ON ASYMPTOMATIC COMPARTMENT 
    gamA = 1/4 
    # RECOVERY RATE ON SYMPTOMATIC NON-TREATMENT COMPARTMENT
    gamY = 1/4 
    # RECOVERY RATE IN HOSPITALIZED COMPARTMENT
    gamH = 1/10.7
    # RATE FROM SYMPTOM ONSET TO HOSPITALIZED
    eta = 0.1695
    # SYMPTOMATIC PROPORTION
    tau = 0.57 
    # EXPOSURE RATE
    sig = 1/2.9
    # PRE-ASYMPTOMATIC COMPARTMENT EXIT RATE
    rhoA = 1/2.3 
    # PRE-SYMLTOMATIC COMPARTMENT EXIT RATE 
    rhoY = 1/2.3 
    # PROPORTION OF PRE-SYMPTOMATIC TRANSMISSION
    P = 0.44
    # RELATIVE INFECTIOUSNESS OF SYMPTOMATIC INDIVIDUALS
    wY = 1.0
    # RELATIVE INFECTIOUSNESS OF INFECTIOUS INDIVIDUALS IN COMPARTMENT I^A
    wA = 0.66
    
    ##############
    # DERIVED SYSTEM PARAMETERS 
    #############

    # INFECTED FATALITY RATIO, AGE SPECIFIC (%)
    IFR = np.array([0.6440/100,6.440/100])
    # SYMPTOMATIC FATALITY RATIO, AGE SPECIFIC (%)
    YFR = np.array([IFR[0]/tau,IFR[1]/tau]) 
    # SYMPTOMATIC CASE HOSPITALIZATION RATE (%)
    YHR = np.array([4.879/100,48.79/100])
    # HOSPITALIZED FATALITY RATIO, AGE SPECIFIC (%)
    HFR = np.array([YFR[0]/YHR[0],YFR[1]/YHR[1]])
    
    # RELATIVE INFECTIOUSNESS OF PRE-SYMTOMATIC INDIVIDUALS
    wP= P/(1-P) /(tau*wY/rhoY + (1-tau)*wA/rhoA) \
        * ((1-tau)*wA/gamA \
        + tau*wY* np.array([YHR[0]/eta+(1-YHR[0])/gamY, \
                          YHR[1]/eta+(1-YHR[1])/gamY])) 
            
    wPY = wP*wY
    wPA = wP*wA 
    
    # RATE OF SYMPTOMATIC INDIVIDUALS GO TO HOSPITAL, AGE-SPECIFIC
    Pi = gamY*np.array([YHR[0]/(eta+(gamY-eta)*YHR[0]),\
                      YHR[1]/(eta+(gamY-eta)*YHR[1])]) # Array with two values
    
    # RATE AT WHICH TERMINAL PATIENTS DIE
    mu = 1/8.1 ###1/mu~8.1

    # TOTAL VENTILATOR CAPACITY IN ALL HOSPITALS
    # #2352 ventilators in Houston (https://www.click2houston.com/health/2020/04/10/texas-medical-center-data-shows-icu-ventilator-capacity-vs-usage-during-coronavirus-outbreak/)
    theta = 3000 

    # TWO TIMES AS MANY PEOPLE NEED VENTILATORS AS THOSE WHO DIE
    rr = 2 

    # DEATH RATE ON HOSPITALIZED INDIVIDUALS, AGE SPECIFIC
    nu = gamH*np.array([HFR[0]/(mu+(gamH-mu)*HFR[0]),\
                      HFR[1]/(mu+(gamH-mu)*HFR[1])])# Array with two values    
          

    ##########
    # COST
    ##########

    a = np.array([[0,2.3,27],[0,2.3,27]]) # Testing costs
    b = np.array([[0,0,40],[0,0,40]]) # Distancing costs
    c = [100,100]     # Cost 5: Opportunity cost for sickness per day (low and high risk)
    d = [500,750]   # Cost 6: Hospitalization cost per day  (low and high risk)
    e = [100000,75000] # Death cost
    f = [2,2] # Cost multiplier of remaining infected
    g = [10000, 10000] # Cost of remaining Susceptible
    nsubDiv = 10 # Number of subdivisions per interval
    
    # MAXIMUM VALUE FOR CONTROL MEASURES
    # [0]: LOW-RISK TESTING RATE 
    # [1]: HIGH-RISK TESTING RATE 
    # [2]: LOW-RISK DISTANCE RATE 
    # [3]: HIGH-RISK DISTANCE RATE 
    uMax = np.array([0.6666,0.6666,.8,.8]) 

    # LOW LEVEL OF CONTROL FOR TESTING/DISTANCING ( AS A FRACTION OF 'uMax')
    frac = [0.3,0.3]


    # Used in enumerative interval algorithm
    cVals = np.array([[0,frac[0]*uMax[0],uMax[0]],
                    [0,frac[0]*uMax[1],uMax[1]],
                    [0,frac[1]*uMax[2],uMax[2]],
                    [0,frac[1]*uMax[3],uMax[3]]])
    
    
 
    params = [wY,wA,wPY,wPA,beta,sig,tau,rhoA,rhoY,
              gamA,gamY,gamH,Pi,eta,nu,mu,theta,rr,Phi]  


    ##############
    # FINAL TIME AND TIME STEP
    ##############

    tf = 180  
    nt = 180


    dictVar={'N':N,
            'NA':NA,
            'xI':xI,
            'beta':beta,
            'Phi':Phi,
            'gamA':gamA,
            'gamY':gamY,
            'gamH':gamH,
            'eta':eta,
            'tau':tau,
            'sig':sig,
            'rhoA':rhoA,
            'rhoY':rhoY,
            'P':P,
            'wY':wY,
            'wA':wA,
            'IFR':IFR,
            'YFR':YFR,
            'YHR':YHR,
            'HFR':HFR,
            'wP':wP,
            'wPY':wPY,
            'wPA':wPA,
            'Pi':Pi,
            'mu':mu,
            'theta':theta,
            'nu':nu,        
            'rr':rr,
            'a':a,
            'b':b,
            'c':c,
            'd':d,
            'e':e,
            'f':f,
            'g':g,
            'tf':tf,
            'nt':nt,
            'nsubDiv':nsubDiv,
            'uMax':uMax,
            'frac':frac,
            'params':params,
            'cVals':cVals}



    return dictVar








