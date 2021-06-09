import numpy as np
import matplotlib.pyplot as plt

from mpl_toolkits.axes_grid1 import make_axes_locatable
from textwrap import wrap


import json


######################################################
### CREATE A LINE GRAPH FOR THE DIFFERENT COST
######################################################
def compartmentPlot(x,costArr,costBestArr,xlabel,ylabel,showPlot = True):
    '''
    Parameter:
        x: n-array for the x-axis
        costArr: (n,9) array
        costBestArr: n-array 
    '''
    
    markers = ['.','^','s','p','P','*','x','D','1','o']
    
    
    plt.yscale('linear')

    nRuns, nCost = costArr.shape

    for ii in range(nCost):
        plt.plot( x, costArr[:,ii])

    plt.plot(x,costBestArr)
    
    plt.legend(['Low Risk Testing',
                'High Risk Testing',
                'Low Risk Distancing',
                'High Risk Distancing',
                'Sickness Opportunity',
                'Hospitalization',
                'Herd',
                'Death',
                'Total Cost'],bbox_to_anchor=(1.04,0.5), loc="center left", borderaxespad=0)   
    
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    
    plt.yscale('log')
    plt.tight_layout()
    
    
    
    if showPlot:
        plt.show()
    

######################################################
### CREATE IMSHOW GRAPHS FOR THE LOW/HIGH TESTING AND DISTANCING COST
######################################################
def uCtrlPlots(uMax,uCtrl,xlabel,ylabel,xmin,xmax,ymin,ymax,showPlot = True):
    

    box = [xmin,xmax,ymin,ymax]
    
    # plt.figure()
    fig2 , ( (ax1,ax2) , (ax3,ax4) )  = plt.subplots(2,2,sharex='col',sharey='row')
    # fig2.set_size_inches(8.0, 6.0)
    
    
    ax1.imshow(uCtrl[:,:,0], 
                origin='lower',
                interpolation='gaussian',
                aspect = 'auto',
                vmin = 0,vmax = uMax[0],
                extent=box)
    #ax1.set_xlabel('Day')
    ax1.set_ylabel(ylabel)
    ax1.set_title('Testing (Low Risk)')
    
    
    
    im2 = ax2.imshow(uCtrl[:,:,1], 
                origin='lower',
                interpolation='gaussian',
                aspect='auto',
                vmin = 0,vmax = uMax[1],
                extent=box)
    #ax2.set_xlabel('Day')
    ax2.set_title('Testing (High Risk)')


    divider = make_axes_locatable(ax2)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    cBar2 = fig2.colorbar(im2,cax = cax,orientation='vertical')
    cBar2.set_label('Control Level')


    ax3.imshow(uCtrl[:,:,2], 
                origin='lower',
                interpolation='gaussian',
                aspect = 'auto',
                vmin = 0,vmax = uMax[2],
                extent=box)
    ax3.set_xlabel('Day')
    ax3.set_ylabel(ylabel)
    ax3.set_title('Social Distance (Low Risk)')
    
    
    im4= ax4.imshow(uCtrl[:,:,3], 
                origin='lower',
                interpolation='gaussian',
                aspect = 'auto',
                vmin = 0,vmax = uMax[3],
                extent=box)
    ax4.set_xlabel('Day')
    ax4.set_title('Social Distance (High Risk)')
    
    divider4 = make_axes_locatable(ax4)
    cax4 = divider4.append_axes('right', size='5%', pad=0.05)
    cBar4 = fig2.colorbar(im4,cax = cax4,orientation='vertical')
    cBar4.set_label('Control Level')

    

    plt.tight_layout()

    if showPlot :
        plt.show()

    return [ax1,ax2,ax3,ax4]


def Heatmaps(filename, dataOnly= False, showPlot = True):


    with open(filename,'r') as f:
        data = json.load(f)

    costArr = np.array(data['costArr'] )
    runtimes = np.array( data['runtime'] )
    xBestArr = np.array( data['xIArr'] )

    uCtrl = np.array( data['uCtrl']['opt'] )

    # get the XY ranges from this data 
    eArr = np.array( data['changes']['eArr'] )[:,0]
    gArr = np.array( data['changes']['gArr'] )[:,0]

    #box = [min(gArr), max(gArr), min(eArr), max(eArr)]
    box = [min(np.log10(gArr)), max(np.log10(gArr)), min(np.log10(eArr)), max(np.log10(eArr))]

    enumCost = np.array( data['totCost']['enum'] )
    mcCost = np.array( data['totCost']['mc'] )
    optCost = np.array( data['totCost']['opt'] )

    uMax = np.array( data['Parameter']['uMax'])
    tf = np.array( data['Parameter']['tf'])



    # FIGURE 1
    totDeath = xBestArr[:,:,-1,8] + xBestArr[:,:,-1,17]
    # FIGURE 2
    implCost = costArr[:,:,0] + costArr[:,:,1] + costArr[:,:,2] + costArr[:,:,3] 
    # FIGURE 3
    nonDeathCost = costArr[:,:,4] + costArr[:,:,5]
    # FIGURE 4
    totNonImmune =  np.sum(xBestArr[:,:,-1,:4] , axis = 2) + np.sum(xBestArr[:,:,-1,9:13],axis = 2) 
    # FIGURE 5
    totRecovered = xBestArr[:,:,-1,7] + xBestArr[:,:,-1,16]
    # FIGURE 6
    totAlive = np.sum(xBestArr[:,:,-1,0:8], axis = 2) + np.sum(xBestArr[:,:,-1,9:17], axis = 2)
    herdImm = totRecovered / totAlive
    # FIGURE 7
    totInfected = np.sum(xBestArr[:,:,-1,5:7], axis = 2) + np.sum(xBestArr[:,:,-1,14:16], axis = 2)
    # FIGURE 8
    totruntime = runtimes[:,:,0] + runtimes[:,:,1] + runtimes[:,:,2]

    plotData = {'totDeath':totDeath,
                'implCost':implCost,
                'nonDeathCost':nonDeathCost,
                'totNonImmune':totNonImmune,
                'totRecovered':totRecovered,
                'hermImm':herdImm,
                'totInfected':totInfected,
                'totRuntime':totruntime,
                'yaxis':eArr,
                'xaxis':gArr,
                'uCtrl':uCtrl,
                'uMax':uMax,
                'tf':tf,
                'enumCost':enumCost,
                'mcCost':mcCost,
                'optCost':optCost}

    if dataOnly:
        return plotData, None

    ####################
    ### Total number of deaths
    ####################

    fig1 = plt.figure(1)
    ax1 = fig1.add_subplot()

    totDeath = np.log10(totDeath)


    im1 = ax1.imshow( totDeath, 
                origin='lower',
                interpolation='gaussian',
                aspect='auto',extent= box )

    #cBar1 = fig1.colorbar(im1, cax = cax1, orientation='vertical' )
    cBar1 = fig1.colorbar(im1, orientation='vertical' )
    cBar1.set_label('$Log_{10}$( Total number of deaths )')

    ax1.set_ylabel('$Log_{10}$( Death cost )')
    ax1.set_xlabel('$Log_{10}$( Non-Immunity cost )')
    ax1.set_title('Total number of deaths',wrap = True)

    plt.tight_layout()


    ####################
    ### Implementation Cost
    ####################

    fig2 = plt.figure(2)
    ax2 = fig2.add_subplot()

    implCost = np.log10(implCost)

    # # Possible zeros
    implCost[ implCost == -np.inf ] = 0

    im2 = ax2.imshow( implCost , 
                origin='lower',
                interpolation='gaussian',
                aspect='auto',extent= box )

    cBar2 = fig2.colorbar(im2, orientation='vertical' )
    cBar2.set_label('$Log_{10}$( Implementation cost ) in \$USD')

    ax2.set_ylabel('$Log_{10}$( Death cost )')
    #plt.xlabel('Susceptible Cost')
    ax2.set_xlabel('$Log_{10}$( Non-Immunity cost )')
    ax2.set_title('Implementation cost',wrap = True)

    plt.tight_layout()

    ####################
    ### Non-Death Impact Cost
    ####################

    fig3 = plt.figure(3)
    ax3 = fig3.add_subplot()

    nonDeathCost = np.log10(nonDeathCost)

    im3 = ax3.imshow( nonDeathCost , 
                origin='lower',
                interpolation='gaussian',
                aspect='auto',extent= box )

    cBar3 = fig3.colorbar(im3, orientation='vertical' )
    cBar3.set_label('$Log_{10}$( Non-death impact cost ) in \$USD')

    ax3.set_ylabel('$Log_{10}$( Death cost )')
    #plt.xlabel('Susceptible Cost')
    ax3.set_xlabel('$Log_{10}$( Non-Immunity cost ) ')
    ax3.set_title('Non-death impact cost',wrap = True)

    plt.tight_layout()

    ####################
    ### Total number of remaining non-immune
    ####################

    fig4 = plt.figure(4)
    ax4 = fig4.add_subplot()

    totNonImmune = np.log10(totNonImmune)

    im4 = ax4.imshow( totNonImmune, 
                origin='lower',
                interpolation='gaussian',
                aspect='auto',extent= box )

    cBar4 = fig4.colorbar(im4, orientation='vertical' )
    cBar4.set_label('$Log_{10}$( Total number of remaining non-immune )')

    ax4.set_ylabel('$Log_{10}$( Death cost )')
    #plt.xlabel('Susceptible Cost')
    ax4.set_xlabel('$Log_{10}$( Non-Immunity cost )')
    ax4.set_title('Total number of remaining non-immune',wrap = True)

    plt.tight_layout()

    ####################
    ### RECOVERED
    ####################

    fig5 = plt.figure(5)
    ax5 = fig5.add_subplot()

    totRecovered = np.log10(totRecovered)

    im5 = ax5.imshow( totRecovered, 
                origin='lower',
                interpolation='gaussian',
                aspect='auto',extent= box )

    cBar5 = fig5.colorbar(im5, orientation='vertical' )
    cBar5.set_label('$Log_{10}$( Total number of recovered )')

    ax5.set_ylabel('$Log_{10}$( Death cost )')
    #plt.xlabel('Susceptible Cost')
    ax5.set_xlabel('$Log_{10}$( Non-Immunity cost )')
    ax5.set_title('Total number of recovered ',wrap = True)


    plt.tight_layout()



    ####################
    ### HERD IMMUNITY
    ####################

    fig6 = plt.figure(6)
    ax6 = fig6.add_subplot()

    im6 = ax6.imshow( herdImm, 
                origin='lower',
                interpolation='gaussian',
                aspect='auto',extent= box )

    cBar6 = fig6.colorbar(im6, orientation='vertical' )
    cBar6.set_label('Level of Herd Immunity')

    ax6.set_ylabel('$Log_{10}$( Death cost )')
    #plt.xlabel('Susceptible Cost')
    ax6.set_xlabel('$Log_{10}$( Non-Immunity cost )')
    ax6.set_title('Herd Immunity ',wrap = True) 


    ####################
    ### Infected
    ####################

    fig7 = plt.figure(7)
    ax7 = fig7.add_subplot()

    im7 = ax7.imshow( np.log10(totInfected), 
                origin='lower',
                interpolation='gaussian',
                aspect='auto',extent= box )

    cBar7 = fig7.colorbar(im7, orientation='vertical' )
    cBar7.set_label('$Log_{10}$( Total Number of remained infected )')

    ax7.set_ylabel('$Log_{10}$( Death cost )')
    #plt.xlabel('Susceptible Cost')
    ax7.set_xlabel('$Log_{10}$( Non-Immunity cost )')
    ax7.set_title('Total number of remaining infected',wrap = True) 





    # ####################
    # ### Execution time
    # ####################
    fig = plt.figure()

    totruntime = runtimes[:,:,0] + runtimes[:,:,1] + runtimes[:,:,2]
    # #totruntime = runtimes[:,:,0] / 60
    # #totruntime = runtimes[:,:,1] / 60
    # #totruntime = runtimes[:,:,2] / 60

    # totruntime = totruntime / 60

    x,y,_ = runtimes.shape
    print('shape:({},{})'.format(x,y))
    print('Runs:{} Runtime:{:.3f} hrs'.format(x*y,np.sum(runtimes)/3600))
    print('----------- Avg Runtime ---------')
    print('enumMean:',np.mean(runtimes[:,:,0] / 60) , ' min' )
    print('mcMean:',np.mean(runtimes[:,:,1] / 60) , ' min')
    print('optMean:',np.mean(runtimes[:,:,2] / 60) , ' min')
    print('Total:',np.mean(totruntime / 60) , ' min')

    im = plt.imshow( totruntime / 60, 
                origin='lower',
                aspect='auto',
                extent= box )

    plt.ylabel('$Log_{10}$( Death cost )')
    plt.xlabel('$Log_{10}$( Non-Immunity cost )')
    plt.title('Execution time',wrap = True)

    cBar = plt.colorbar(im)
    cBar.set_label('Exec. runtime (min)')

    plt.tight_layout()

    if showPlot:
        plt.show()

    return plotData, [ [fig1,ax1],[fig2,ax2],[fig3,ax3],[fig4,ax4],[fig5,ax5],[fig6,ax6],[fig7,ax7] ]


def pareto_plot(filename,deathArr,recArr, dataOnly = False ):


    # GET THE FIGURES CREATED IN THIS METHOD
    plotData, axesList = Heatmaps(filename, dataOnly = dataOnly, showPlot = False)

    # minimize implementation cost
    with open(filename,'r') as f:
        data = json.load(f)

    uMax = np.array( data['Parameter']['uMax'])
    tf = np.array( data['Parameter']['tf'])

    costArr = np.array(data['costArr'] )
    uCtrl = np.array( data['uCtrl']['opt'] )

    # # get the XY ranges from this data 
    eArr = np.array( data['changes']['eArr'] )[:,0]
    gArr = np.array( data['changes']['gArr'] )[:,0]

    # FIGURE 1
    totDeath = plotData['totDeath']
    # FIGURE 2
    implCost = plotData['implCost']
    # FIGURE 3
    nonDeathCost = plotData['nonDeathCost']
    # FIGURE 4
    totNonImmune =  plotData['totNonImmune']
    # FIGURE 5
    totRecovered = plotData['totRecovered']
    # FIGURE 6
    herdImm = plotData['hermImm']

    enumCost = np.array( data['totCost']['enum'] )
    mcCost = np.array( data['totCost']['mc'] )
    optCost = np.array( data['totCost']['opt'] )

    # LOG10 VERSIONS
    totDeath_Log10 = np.log10( np.maximum( totDeath , 0.09 ) )
    implCost_Log10 = np.log10( np.maximum( implCost , 0.09 ) )
    nonDeathCost_Log10 = np.log10( np.maximum( nonDeathCost , 0.09))
    totNonImmune_Log10 = np.log10( np.maximum( totNonImmune , 0.09))
    totRecovered_Log10 = np.log10( np.maximum( totRecovered , 0.09))


    ###################
    ### PARETO CODE
    ###################

    valueArr = np.zeros( ( len(deathArr) , 15 ) ) 

    vBool = np.zeros( len(deathArr)  , dtype = bool)
    dataBool = np.zeros( (len(eArr),len(gArr) ) , dtype = bool)

    #gridSize = totDeath.shape[0] * totDeath.shape[1]

#    for ii,rr in enumerate(recArr):

#        bool1 = totRecovered >= rr

    for jj,dd in enumerate(deathArr):
     
        bool2 = totDeath <= dd

        # SELECT STRATEGIES THAT HAVE (RECOVERED >= RR) AND (DEATHS <= DD) 
        indicator = bool2

        #  CHECK THAT ATLEAST ONE STRATEGY WAS SELECTED AND NOT ALL STRATEGIES GET SELECTED
        if np.sum(indicator) > 0 :#and np.sum(indicator) < gridSize :

            # FROM THE SELECTED STRATEGY, SELECT THOSE WITH THE MAX NUMBER OF RECOVERED
            bool5 = ( totRecovered == np.max(totRecovered[indicator]) )
            indicator = np.logical_and(indicator,bool5)

            tmp = 1.0 * implCost
            notIt = np.invert(indicator) # SELECT STRATEGIES THAT DID NOT MEET CRITERIA
            tmp[notIt] = np.inf # SET THE NON-SELECTED STRATEGIES' IMPLEMENTATION COST TO INF

            #tmp[dataBool] = np.inf # SET THE IMPLEMENTATION COST OF STRATEGIES ALREADY RECORDED TO INF

            # SELECT STRATEGIES THAT MINIMIZE THE IMPLEMENTATION COST
            bool6 = ( tmp == np.min(tmp) )
            indicator = np.logical_and(indicator,bool6) 

            # # FROM THE SELECTED STRATEGY, SELECT THOSE WITH THE MIN NUMBER OF DEATHS
            bool4 = ( totDeath == np.min(totDeath[indicator]) )
            indicator = np.logical_and(indicator,bool4)

            # SELECT INDECES WHERE CRITERIA IS MET ( IF THERE'S MANY, CHOOSE THE FIRST ONE )
            ix = np.argwhere(indicator)[0]

            # IF THIS STRATEGY HAS NOT BEEN SELECTED YET, THEN CONTINUE
            if not(dataBool[ ix[0], ix[1] ]) :

                vBool[jj] = True
                dataBool[ ix[0], ix[1] ] = True

                ## ACTUAL VALUES FROM DATA
                cRecov = totRecovered[ ix[0],ix[1] ]
                cDeath = totDeath[ ix[0], ix[1] ]
                cNonImmune = totNonImmune[ ix[0], ix[1] ]
                cminCost = implCost[ ix[0], ix[1] ]
                cNonDeathCost = nonDeathCost[ ix[0], ix[1] ]
                cImm = herdImm[ix[0],ix[1]]

                valueArr[jj][0] = 0
                valueArr[jj][1] = dd
                valueArr[jj][2] = cDeath 
                valueArr[jj][3] = cRecov 
                valueArr[jj][4] = cNonImmune 
                valueArr[jj][5] = cminCost 
                valueArr[jj][6] = cNonDeathCost 
                valueArr[jj][7] = cImm 

                valueArr[jj][8] = gArr[ix[1]] 
                valueArr[jj][9] = eArr[ix[0]] 

                valueArr[jj][10] = ix[0]
                valueArr[jj][11] = ix[1] 


                valueArr[jj][12] = enumCost[ ix[0], ix[1] ]
                valueArr[jj][13] = mcCost[ ix[0], ix[1] ]
                valueArr[jj][14] = optCost[ ix[0], ix[1] ]


    # SELECT THE ELEMENTS THAT ACTUALLY HAVE DATA IN THEM
    valueVec = valueArr[vBool]


    if dataOnly:
        return valueVec, plotData



    # PLACE THE MARKERS ON THE HEAT MAPS
    for fig, ax in axesList:

        cb = ax.images[0].colorbar 
        label = cb._label
        cb.remove()


        scc = ax.scatter( np.log10(valueVec[:,8] ),np.log10(valueVec[:,9]), s = 100, c = np.log10(valueVec[:,3])  , 
                    vmin = np.min(totRecovered_Log10) , vmax = np.max(totRecovered_Log10) , marker = 'o',edgecolor = 'k')

        cBar = fig.colorbar(scc, orientation = 'vertical')
        cBar.set_label('$Log_{10}$( Tot. number of recovered )')

        cBar1 = fig.colorbar(ax.images[0], orientation = 'vertical')
        cBar1.set_label(label)

        fig.tight_layout()


    ##############
    ### Implementation cost vs Number of Deaths
    #############

    fig1 = plt.figure()
    ax1 = fig1.add_subplot()

    div1 = make_axes_locatable(ax1)
    cax1 = div1.append_axes("right", size='5%', pad=0.05)    


    sc = ax1.scatter( valueVec[:,2], valueVec[:,5], c = np.log10(valueVec[:,3])  , 
                    vmin = np.min(totRecovered_Log10) , vmax = np.max(totRecovered_Log10)  ,  
                    marker = 'o')

    cBar1 = plt.colorbar(sc, cax = cax1 , orientation = 'vertical')
    cBar1.set_label('$Log_{10}$( Total number of recovered )')
    ax1.set_xlabel('Tot. number of deaths')
    ax1.set_ylabel('Implementation Cost in \$USD')

    plt.tight_layout()

    ##############
    ### Numbers of death vs Control measures
    #############

    uCtrlMx = np.zeros( ( len(valueVec[:,2]) , tf, 4 ) )
    compCost = np.zeros( ( len(valueVec[:,2]) , 8 ) )

    for ii in range( len(valueVec) ):

        ee,gg = valueVec[ii,10:12]
        uCtrlMx[ii,:,:] = uCtrl[ int(ee), int(gg) , : ,:]

        compCost[ii,:] = costArr[ int(ee), int(gg) , :]

    idxDeathSort = np.argsort(valueVec[:,2])

    uCtrlPlots(uMax,uCtrlMx[idxDeathSort],'Days','Tot. number of deaths',0,tf,np.min(valueVec[:,2]),np.max(valueVec[:,2]) ,showPlot=False )
    plt.tight_layout()

    ##############
    ### Numbers of death vs Tot. Cost
    #############

    plt.figure()

    deathSort = valueVec[:,2][idxDeathSort]
    enumSort = valueVec[:,12][idxDeathSort]
    mcSort = valueVec[:,13][idxDeathSort]
    localSort = valueVec[:,14][idxDeathSort]

    if np.sum(valueVec[:,12][idxDeathSort]) > 0 :
        plt.plot( deathSort , enumSort / localSort , label ='Enum' , linestyle = '-',linewidth = 2 , color = 'b')

    if np.sum(valueVec[:,13][idxDeathSort]) > 0 :
        plt.plot( deathSort , mcSort / localSort , label = 'Mc' , linestyle ='-', linewidth = 2, color = 'g')

    plt.plot( deathSort , localSort / localSort , label ='Local Opt', linestyle =':',linewidth = 2, color = '#FFA500')

    #plt.yscale('log')
    plt.xlabel('Total number of deaths')
    plt.ylabel('Total Cost in \$USD proportion to local optimization')
    plt.legend()
    plt.tight_layout()


    ##############
    ### Numbers of death vs Compartment Cost
    #############

    plt.figure()

    compartmentPlot(valueVec[:,2][idxDeathSort], compCost[idxDeathSort], valueVec[:,14][idxDeathSort], 'Total number of deaths', 'Total Cost in \$USD',showPlot=False)
    plt.tight_layout()

    plt.show()


def singleRun_plots(filename,rowIdx = 5,colIdx = 2):

    with open(filename,'r') as f:
        data = json.load(f)

    costArr = np.array(data['costArr'] )
    runtimes = np.array( data['runtime'] )
    xBestArr = np.array( data['xIArr'] )

    uCtrl = np.array( data['uCtrl']['opt'] )

    # get the XY ranges from this data 
    eArr = np.array( data['changes']['eArr'] )[:,0]
    gArr = np.array( data['changes']['gArr'] )[:,0]

    #box = [min(gArr), max(gArr), min(eArr), max(eArr)]
    box = [min(np.log10(gArr)), max(np.log10(gArr)), min(np.log10(eArr)), max(np.log10(eArr))]

    enumCost = np.array( data['totCost']['enum'] )
    mcCost = np.array( data['totCost']['mc'] )
    optCost = np.array( data['totCost']['opt'] )


    print('Final Population')
    for pop in xBestArr[rowIdx,colIdx,-1,:] :
        print(pop)

    plt.figure()
    size , _ = uCtrl[rowIdx,colIdx].shape

    plt.plot( range(size) , uCtrl[rowIdx,colIdx,:,0] , label = 'Low-Risk Testing')
    plt.plot( range(size) , uCtrl[rowIdx,colIdx,:,1] , label = 'High-Risk Testing')
    plt.plot( range(size) , uCtrl[rowIdx,colIdx,:,2] , label = 'Low-Risk Distancing')
    plt.plot( range(size) , uCtrl[rowIdx,colIdx,:,3] , label = 'High-Risk Distancing')

    plt.legend()
    plt.title('Local Opt. Control Measures')
    plt.xlabel('Days')
    plt.ylabel('Control rate')

    plt.show()




if __name__ == '__main__':


    # MODIFIED LOCAL OPT
    filename = 'nonImm_Death_partialEnum_partialMC.json'
    #filename = 'nonImm_Death_fullEnum_noMC.json'
    #filename = 'nonImm_Death_noEnum_fullMC.json'
    

    ###########
    # PLOTS
    ###########

    #Heatmaps(filename, dataOnly=True)

    deathArr = np.arange(0,100E3,10E3) 
    recArr = 10 ** np.linspace(2,6,10)#np.array([2,3,4,5,6]) 
    pareto_plot(filename,deathArr,recArr)


    #singleRun_plots(filename,rowIdx = 5,colIdx = 2)


















