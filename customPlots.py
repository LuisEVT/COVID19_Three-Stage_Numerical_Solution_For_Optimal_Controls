import numpy as np
import matplotlib.pyplot as plt
from tbPlot import pareto_plot

from mpl_toolkits.axes_grid1 import ImageGrid
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.ticker as ticker

##########
# THE GRAPHS IN THIS FILE ARE CUSTOM MADE FOR THE PUBLISH PAPER. 
##########


# STATE THE THREE FILES
file_PEPM = 'nonImm_Death_partialEnum_partialMC.json'
file_FENM = 'nonImm_Death_fullEnum_noMC.json'
file_NEFM = 'nonImm_Death_noEnum_FullMC.json'

# USED FOR PARETO
#nDeathArr = np.arange(1E3,90E3,10E3) 
#nRecovArr = 10 ** np.array([3,4,5,6]) 

nDeathArr = np.arange(0,85E3,3E3) 
nRecovArr = 10 ** np.array([2,3,4,5,6,7]) # I'll delete this later


#####################################################
### READ FILES AND RETRIEVE THE INFORMATION NEEDED
#####################################################


# GET THE INFORMATION FOR THE GRAPHS
paretoPts_PEPM, plotData_PEPM = pareto_plot(file_PEPM, nDeathArr, nRecovArr, dataOnly=True)
paretoPts_FENM, plotData_FENM = pareto_plot(file_FENM, nDeathArr, nRecovArr, dataOnly=True)
paretoPts_NEFM, plotData_NEFM = pareto_plot(file_NEFM, nDeathArr, nRecovArr, dataOnly=True)

eArr = plotData_PEPM['yaxis']
gArr = plotData_PEPM['xaxis']
box = [ min(np.log10(gArr)), max(np.log10(gArr)), min(np.log10(eArr)), max(np.log10(eArr)) ]

tf = plotData_PEPM['tf']
uMax = plotData_PEPM['uMax']

# TOTAL NUMBER OF DEATHS
totDeath_PEPM = np.log10( np.maximum( plotData_PEPM['totDeath'] , 0.9 ) )
totDeath_FENM = np.log10( np.maximum( plotData_FENM['totDeath'] , 0.9 ) )
totDeath_NEFM = np.log10( np.maximum( plotData_NEFM['totDeath'] , 0.9 ) )

deathVmin = np.min( [ np.min(totDeath_PEPM) , np.min(totDeath_FENM) , np.min(totDeath_NEFM) ]  )
deathVmax = np.max( [ np.max(totDeath_PEPM) , np.max(totDeath_FENM) , np.max(totDeath_NEFM) ]  )

# GET THE RECOVERED VALUES
totRecovered_PEPM = np.log10( np.maximum( plotData_PEPM['totRecovered'] , 0.9 ) )
totRecovered_FENM = np.log10( np.maximum( plotData_FENM['totRecovered'] , 0.9 ) )
totRecovered_NEFM = np.log10( np.maximum( plotData_NEFM['totRecovered'] , 0.9 ) )

recVmin = np.min( [ np.min(totRecovered_PEPM) , np.min(totRecovered_FENM) , np.min(totRecovered_NEFM) ]  )
recVmax = np.max( [ np.max(totRecovered_PEPM) , np.max(totRecovered_FENM) , np.max(totRecovered_NEFM) ]  )

# GET IMPLEMENTATION COST
implementCost_PEPM = np.log10( np.maximum( plotData_PEPM['implCost'] , 0.9 ) )
implementCost_FENM = np.log10( np.maximum( plotData_FENM['implCost'] , 0.9 ) )
implementCost_NEFM = np.log10( np.maximum( plotData_NEFM['implCost'] , 0.9 ) )

implCostVmin = np.min( [ np.min(implementCost_PEPM) , np.min(implementCost_FENM) , np.min(implementCost_NEFM) ]  )
implCostVmax = np.max( [ np.max(implementCost_PEPM) , np.max(implementCost_FENM) , np.max(implementCost_NEFM) ]  )

# TOTAL NUMBER OF NON-IMMUNE 
totNonImmune_PEPM = np.log10( np.maximum( plotData_PEPM['totNonImmune'] , 0.9 ) )
totNonImmune_FENM = np.log10( np.maximum( plotData_FENM['totNonImmune'] , 0.9 ) )
totNonImmune_NEFM = np.log10( np.maximum( plotData_NEFM['totNonImmune'] , 0.9 ) )

nonImmuneVmin = np.min( [ np.min(totNonImmune_PEPM) , np.min(totNonImmune_FENM) , np.min(totNonImmune_NEFM) ]  )
nonImmuneVmax = np.max( [ np.max(totNonImmune_PEPM) , np.max(totNonImmune_FENM) , np.max(totNonImmune_NEFM) ]  )

finalCost_PEPM = plotData_PEPM['optCost']
finalCost_FENM = plotData_FENM['optCost']
finalCost_NEFM = plotData_NEFM['optCost']

#####################################################################
####################################################################
####################################################################

########
# FIGURE 1: TOTAL NUMBER OF DEATHS
########

fig1, (ax11,ax12,ax13) = plt.subplots(nrows=1, ncols=3 ,sharey = True, sharex = True,figsize=(12, 4), gridspec_kw={'width_ratios': [1,1,1.5]})

im11 = ax11.imshow( totDeath_PEPM, origin='lower', interpolation='gaussian', aspect='auto', vmin= deathVmin, vmax= deathVmax, extent = box)
im12 = ax12.imshow( totDeath_FENM, origin='lower', interpolation='gaussian', aspect='auto', vmin= deathVmin, vmax= deathVmax, extent = box)
im13 = ax13.imshow( totDeath_NEFM, origin='lower', interpolation='gaussian', aspect='auto', vmin= deathVmin, vmax= deathVmax, extent = box)


scc = ax11.scatter( np.log10(paretoPts_PEPM[:,8] ), np.log10(paretoPts_PEPM[:,9]), s = 100, c = np.log10(paretoPts_PEPM[:,3])  , 
            vmin = recVmin , vmax = recVmax , marker = 'o',edgecolor = 'k')

scc = ax12.scatter( np.log10(paretoPts_FENM[:,8] ), np.log10(paretoPts_FENM[:,9]), s = 100, c = np.log10(paretoPts_FENM[:,3])  , 
            vmin = recVmin , vmax = recVmax , marker = 'o',edgecolor = 'k')

scc = ax13.scatter( np.log10(paretoPts_NEFM[:,8] ), np.log10(paretoPts_NEFM[:,9]), s = 100, c = np.log10(paretoPts_NEFM[:,3])  , 
            vmin = recVmin , vmax = recVmax , marker = 'o',edgecolor = 'k')            

cBar = fig1.colorbar(scc, orientation = 'vertical')
cBar.set_label('$Log_{10}$( Total number of recovered )')

cBar1 = fig1.colorbar(im13, orientation='vertical')
cBar1.set_label('$Log_{10}$( Total number of deaths )')

ax11.set_title('PEPM')
ax12.set_title('FENM')
ax13.set_title('NEFM')

ax11.set_ylabel('$Log_{10}$( Death cost )')
ax11.set_xlabel('$Log_{10}$( Non-Immunity cost )')
ax12.set_xlabel('$Log_{10}$( Non-Immunity cost )')
ax13.set_xlabel('$Log_{10}$( Non-Immunity cost )')

fig1.suptitle('Total number of deaths', wrap = True, fontsize=16)

plt.tight_layout()
plt.subplots_adjust(top = 0.83,wspace=0.02)


####################
# Figure 2: Implementation Cost
####################

fig2, (ax21,ax22,ax23) = plt.subplots(nrows=1, ncols=3 ,sharey = True, sharex = True,figsize=(12, 4), gridspec_kw={'width_ratios': [1,1,1.5]})


im21 = ax21.imshow( implementCost_PEPM, origin='lower', interpolation='gaussian', aspect='auto', vmin= implCostVmin, vmax= implCostVmax, extent = box)
im22 = ax22.imshow( implementCost_FENM, origin='lower', interpolation='gaussian', aspect='auto', vmin= implCostVmin, vmax= implCostVmax, extent = box)
im23 = ax23.imshow( implementCost_NEFM, origin='lower', interpolation='gaussian', aspect='auto', vmin= implCostVmin, vmax= implCostVmax, extent = box)


scc = ax21.scatter( np.log10(paretoPts_PEPM[:,8] ), np.log10(paretoPts_PEPM[:,9]), s = 100, c = np.log10(paretoPts_PEPM[:,3])  , 
            vmin = recVmin , vmax = recVmax , marker = 'o',edgecolor = 'k')

scc = ax22.scatter( np.log10(paretoPts_FENM[:,8] ), np.log10(paretoPts_FENM[:,9]), s = 100, c = np.log10(paretoPts_FENM[:,3])  , 
            vmin = recVmin , vmax = recVmax , marker = 'o',edgecolor = 'k')

scc = ax23.scatter( np.log10(paretoPts_NEFM[:,8] ), np.log10(paretoPts_NEFM[:,9]), s = 100, c = np.log10(paretoPts_NEFM[:,3])  , 
            vmin = recVmin , vmax = recVmax , marker = 'o',edgecolor = 'k')            

cBar = fig2.colorbar(scc, orientation = 'vertical')
cBar.set_label('$Log_{10}$( Total number of recovered )')



cBar2 = fig2.colorbar(im23, orientation='vertical')
cBar2.set_label('$Log_{10}$( Total Implementation cost )')

ax21.set_title('PEPM')
ax22.set_title('FENM')
ax23.set_title('NEFM')

ax21.set_ylabel('$Log_{10}$( Death cost )')
ax21.set_xlabel('$Log_{10}$( Non-Immunity cost )')
ax22.set_xlabel('$Log_{10}$( Non-Immunity cost )')
ax23.set_xlabel('$Log_{10}$( Non-Immunity cost )')

fig2.suptitle('Total Implementation cost', wrap = True, fontsize=16)

plt.tight_layout()
plt.subplots_adjust(top = 0.83,wspace=0.02)


####################
# Figure 3: Total number of remaining non-immune
####################

fig3, (ax31,ax32,ax33) = plt.subplots(nrows=1, ncols=3 ,sharey = True, sharex = True,figsize=(12, 4), gridspec_kw={'width_ratios': [1,1,1.5]})

im31 = ax31.imshow( totNonImmune_PEPM, origin='lower', interpolation='gaussian', aspect='auto', vmin= nonImmuneVmin, vmax= nonImmuneVmax, extent = box)
im32 = ax32.imshow( totNonImmune_FENM, origin='lower', interpolation='gaussian', aspect='auto', vmin= nonImmuneVmin, vmax= nonImmuneVmax, extent = box)
im33 = ax33.imshow( totNonImmune_NEFM, origin='lower', interpolation='gaussian', aspect='auto', vmin= nonImmuneVmin, vmax= nonImmuneVmax, extent = box)


scc = ax31.scatter( np.log10(paretoPts_PEPM[:,8] ), np.log10(paretoPts_PEPM[:,9]), s = 100, c = np.log10(paretoPts_PEPM[:,3])  , 
            vmin = recVmin , vmax = recVmax , marker = 'o',edgecolor = 'k')

scc = ax32.scatter( np.log10(paretoPts_FENM[:,8] ), np.log10(paretoPts_FENM[:,9]), s = 100, c = np.log10(paretoPts_FENM[:,3])  , 
            vmin = recVmin , vmax = recVmax , marker = 'o',edgecolor = 'k')

scc = ax33.scatter( np.log10(paretoPts_NEFM[:,8] ), np.log10(paretoPts_NEFM[:,9]), s = 100, c = np.log10(paretoPts_NEFM[:,3])  , 
            vmin = recVmin , vmax = recVmax , marker = 'o',edgecolor = 'k')            

cBar = fig3.colorbar(scc, orientation = 'vertical')
cBar.set_label('$Log_{10}$( Total number of recovered )')


cBar3 = fig3.colorbar(im33, orientation='vertical')
cBar3.set_label('$Log_{10}$( Total number of remaining non-immune )')

ax31.set_title('PEPM')
ax32.set_title('FENM')
ax33.set_title('NEFM')

ax31.set_ylabel('$Log_{10}$( Death cost )')
ax31.set_xlabel('$Log_{10}$( Non-Immunity cost )')
ax32.set_xlabel('$Log_{10}$( Non-Immunity cost )')
ax33.set_xlabel('$Log_{10}$( Non-Immunity cost )')

fig3.suptitle('Total number of remaining non-immune', wrap = True, fontsize=16)

plt.tight_layout()
plt.subplots_adjust(top = 0.83,wspace=0.02)


####################
# Figure 4: Total number of Recovered
####################

fig4, (ax41,ax42,ax43) = plt.subplots(nrows=1, ncols=3 ,sharey = True, sharex = True,figsize=(12, 4), gridspec_kw={'width_ratios': [1,1,1.5]})


im41 = ax41.imshow( totRecovered_PEPM, origin='lower', interpolation='gaussian', aspect='auto', vmin= recVmin, vmax= recVmax, extent = box)
im42 = ax42.imshow( totRecovered_FENM, origin='lower', interpolation='gaussian', aspect='auto', vmin= recVmin, vmax= recVmax, extent = box)
im43 = ax43.imshow( totRecovered_NEFM, origin='lower', interpolation='gaussian', aspect='auto', vmin= recVmin, vmax= recVmax, extent = box)

scc = ax41.scatter( np.log10(paretoPts_PEPM[:,8] ), np.log10(paretoPts_PEPM[:,9]), s = 100, c = np.log10(paretoPts_PEPM[:,3])  , 
            vmin = recVmin , vmax = recVmax , marker = 'o',edgecolor = 'k')

scc = ax42.scatter( np.log10(paretoPts_FENM[:,8] ), np.log10(paretoPts_FENM[:,9]), s = 100, c = np.log10(paretoPts_FENM[:,3])  , 
            vmin = recVmin , vmax = recVmax , marker = 'o',edgecolor = 'k')

scc = ax43.scatter( np.log10(paretoPts_NEFM[:,8] ), np.log10(paretoPts_NEFM[:,9]), s = 100, c = np.log10(paretoPts_NEFM[:,3])  , 
            vmin = recVmin , vmax = recVmax , marker = 'o',edgecolor = 'k')            

cBar = fig4.colorbar(scc, orientation = 'vertical')
cBar.set_label('$Log_{10}$( Total number of recovered )')

cBar4 = fig4.colorbar(im43, orientation='vertical')
cBar4.set_label('$Log_{10}$( Total number of recovered )')

ax41.set_title('PEPM')
ax42.set_title('FENM')
ax43.set_title('NEFM')

ax41.set_ylabel('$Log_{10}$( Death cost )')
ax41.set_xlabel('$Log_{10}$( Non-Immunity cost )')
ax42.set_xlabel('$Log_{10}$( Non-Immunity cost )')
ax43.set_xlabel('$Log_{10}$( Non-Immunity cost )')

fig4.suptitle('Total number of recovered', wrap = True, fontsize=16)

plt.tight_layout()
plt.subplots_adjust(top = 0.83,wspace=0.02)



####################
# Figure 5: Execution Time
####################

fig5, (ax51,ax52,ax53) = plt.subplots(nrows=1, ncols=3 ,sharey = True, sharex = True,figsize=(12, 4), gridspec_kw={'width_ratios': [1,1,1.2]})

totRuntime_PEPM = plotData_PEPM['totRuntime'] / 60
totRuntime_FENM = plotData_FENM['totRuntime'] / 60
totRuntime_NEFM = plotData_NEFM['totRuntime'] / 60

fig5_vmin = 0
fig5_vmax = 120

im51 = ax51.imshow( totRuntime_PEPM, origin='lower', aspect='auto', vmin= fig5_vmin, vmax= fig5_vmax, extent = box)
im52 = ax52.imshow( totRuntime_FENM, origin='lower', aspect='auto', vmin= fig5_vmin, vmax= fig5_vmax, extent = box)
im53 = ax53.imshow( totRuntime_NEFM, origin='lower', aspect='auto', vmin= fig5_vmin, vmax= fig5_vmax, extent = box)

cBar5 = fig5.colorbar(im53, orientation='vertical')
cBar5.set_label(' Execution time (min) ')

ax51.set_title('PEPM')
ax52.set_title('FENM')
ax53.set_title('NEFM')

ax51.set_ylabel('$Log_{10}$( Death cost )')
ax51.set_xlabel('$Log_{10}$( Non-Immunity cost )')
ax52.set_xlabel('$Log_{10}$( Non-Immunity cost )')
ax53.set_xlabel('$Log_{10}$( Non-Immunity cost )')

fig5.suptitle('Execution time', wrap = True, fontsize=16)

plt.tight_layout()
plt.subplots_adjust(top = 0.83,wspace=0.02)


####################
# Figure 6.1: Implementation cost vs Number of Deaths ( W/ Recovery)
####################

fig6_1, ax6_1 = plt.subplots(figsize=(6, 4))

xSort_PEPM = np.argsort(paretoPts_PEPM[:,2])
xSort_FENM = np.argsort(paretoPts_FENM[:,2])
xSort_NEFM = np.argsort(paretoPts_NEFM[:,2])

lp = ax6_1.plot(paretoPts_PEPM[:,2][xSort_PEPM], paretoPts_PEPM[:,5][xSort_PEPM], linewidth = 2, label = 'PEPM',zorder=1)    
lp = ax6_1.plot(paretoPts_FENM[:,2][xSort_FENM], paretoPts_FENM[:,5][xSort_FENM], linewidth = 2, label = 'FENM',zorder=1)
lp = ax6_1.plot(paretoPts_NEFM[:,2][xSort_NEFM], paretoPts_NEFM[:,5][xSort_NEFM], linewidth = 2, label = 'NEFM',zorder=1) 

sc = ax6_1.scatter( paretoPts_PEPM[:,2], paretoPts_PEPM[:,5], c = np.log10(paretoPts_PEPM[:,3])  , 
                vmin= recVmin, vmax= recVmax,  marker = 'o',zorder=2 )

sc = ax6_1.scatter( paretoPts_FENM[:,2], paretoPts_FENM[:,5], c = np.log10(paretoPts_FENM[:,3])  , 
                vmin= recVmin, vmax= recVmax,  marker = 'o',zorder=2)           

sc = ax6_1.scatter( paretoPts_NEFM[:,2], paretoPts_NEFM[:,5], c = np.log10(paretoPts_NEFM[:,3])  , 
                vmin= recVmin, vmax= recVmax,  marker = 'o',zorder=2)

divider = make_axes_locatable(ax6_1)
cax = divider.append_axes('right', size='5%', pad=0.05)
cBar1 = plt.colorbar(sc, orientation = 'vertical', cax = cax)
cBar1.set_label('$Log_{10}$(Total number of recovered)\n')
cBar1.ax.set_yticklabels(["{:.2}".format(i) for i in cBar1.get_ticks()]) # set ticks of your format
 

ax6_1.legend()
ax6_1.set_xlabel('Total number of deaths')
ax6_1.set_ylabel('Implementation Cost in \$USD')
plt.tight_layout()

####################
# Figure 6.2: Non-Immune vs Death Pareto Graph ( W/ Cost )
####################

fig6_2, ax6_2 = plt.subplots(figsize=(6, 4))

xSort_PEPM = np.argsort(paretoPts_PEPM[:,2])
xSort_FENM = np.argsort(paretoPts_FENM[:,2])
xSort_NEFM = np.argsort(paretoPts_NEFM[:,2])

paretoImplCostVmin = np.log10(np.min([np.min(paretoPts_PEPM[:,5]),np.min(paretoPts_FENM[:,5]),np.min(paretoPts_NEFM[:,5])] ) )
paretoImplCostVmax = np.log10(np.max([np.max(paretoPts_PEPM[:,5]),np.max(paretoPts_FENM[:,5]),np.max(paretoPts_NEFM[:,5])] ) )


lp = ax6_2.plot(paretoPts_PEPM[:,2][xSort_PEPM], np.log10(paretoPts_PEPM[:,4][xSort_PEPM]), linewidth = 2, label = 'PEPM',zorder=1)    
lp = ax6_2.plot(paretoPts_FENM[:,2][xSort_FENM], np.log10(paretoPts_FENM[:,4][xSort_FENM]), linewidth = 2, label = 'FENM',zorder=1)
lp = ax6_2.plot(paretoPts_NEFM[:,2][xSort_NEFM], np.log10(paretoPts_NEFM[:,4][xSort_NEFM]), linewidth = 2, label = 'NEFM',zorder=1) 

sc = ax6_2.scatter( paretoPts_PEPM[:,2], np.log10(paretoPts_PEPM[:,4]), c = np.log10(paretoPts_PEPM[:,5])  , 
                vmin= paretoImplCostVmin , vmax= paretoImplCostVmax,  marker = 'o',zorder=2)    

sc = ax6_2.scatter( paretoPts_FENM[:,2], np.log10(paretoPts_FENM[:,4]), c = np.log10(paretoPts_FENM[:,5])  , 
                vmin= paretoImplCostVmin, vmax= paretoImplCostVmax,  marker = 'o',zorder=2)                   

sc = ax6_2.scatter( paretoPts_NEFM[:,2], np.log10(paretoPts_NEFM[:,4]), c = np.log10(paretoPts_NEFM[:,5])  , 
                vmin= paretoImplCostVmin, vmax= paretoImplCostVmax,  marker = 'o',zorder=2)    

divider = make_axes_locatable(ax6_2)
cax = divider.append_axes('right', size='5%', pad=0.05)
cBar1 = plt.colorbar(sc, orientation = 'vertical', cax = cax)
cBar1.set_label('$Log_{10}$(Implementation Cost) in \$USD')
cBar1.ax.set_yticklabels(["{:.2}".format(i) for i in cBar1.get_ticks()]) # set ticks of your format

ax6_2.legend()
ax6_2.set_xlabel('Total number of deaths')
ax6_2.set_ylabel('$Log_{10}$(Total remaining non-immune)')
plt.tight_layout()

####################
# Figure 6.3: Implementation Cost vs Non-Immune Pareto Graph (W/ Deaths)
####################

fig6_3, ax6_3 = plt.subplots(figsize=(6, 4))

xSort_PEPM = np.argsort(paretoPts_PEPM[:,4])
xSort_FENM = np.argsort(paretoPts_FENM[:,4])
xSort_NEFM = np.argsort(paretoPts_NEFM[:,4])


lp = ax6_3.plot(np.log10(paretoPts_PEPM[:,4][xSort_PEPM]), paretoPts_PEPM[:,5][xSort_PEPM], linewidth = 2, label = 'PEPM',zorder=1)    
lp = ax6_3.plot(np.log10(paretoPts_FENM[:,4][xSort_FENM]), paretoPts_FENM[:,5][xSort_FENM], linewidth = 2, label = 'FENM',zorder=1)
lp = ax6_3.plot(np.log10(paretoPts_NEFM[:,4][xSort_NEFM]), paretoPts_NEFM[:,5][xSort_NEFM], linewidth = 2, label = 'NEFM',zorder=1) 

sc = ax6_3.scatter( np.log10(paretoPts_PEPM[:,4]), paretoPts_PEPM[:,5], c = np.log10(paretoPts_PEPM[:,2])  , 
                vmin= deathVmin, vmax= deathVmax,  marker = 'o',zorder=2)    

sc = ax6_3.scatter( np.log10(paretoPts_FENM[:,4]), paretoPts_FENM[:,5], c = np.log10(paretoPts_FENM[:,2])  , 
                vmin= deathVmin, vmax= deathVmax,  marker = 'o',zorder=2)                 

sc = ax6_3.scatter( np.log10(paretoPts_NEFM[:,4]), paretoPts_NEFM[:,5], c = np.log10(paretoPts_NEFM[:,2])  , 
                vmin= deathVmin, vmax= deathVmax,  marker = 'o',zorder=2)    


divider = make_axes_locatable(ax6_3)
cax = divider.append_axes('right', size='5%', pad=0.05)
cBar1 = plt.colorbar(sc, orientation = 'vertical', cax = cax)
cBar1.set_label('$Log_{10}$(Total number of deaths)')
cBar1.ax.set_yticklabels(["{:.2}".format(i) for i in cBar1.get_ticks()]) # set ticks of your format


ax6_3.legend()


#cBar1 = plt.colorbar(sc, orientation = 'vertical')
#cBar1.set_label('$Log_{10}$( Total number of recovered )')






ax6_3.set_xlabel('$Log_{10}$(Total remaining non-immune)')
ax6_3.set_ylabel('Implementation Cost in \$USD')

# ax6_1.set_title('PEPM')
# ax6_2.set_title('FENM')
# ax6_3.set_title('NEFM')

# fig6.suptitle('Pareto graphs' ,fontsize=16)

plt.tight_layout()
#plt.subplots_adjust(top = 0.85)


####################
# Figure 6X: Implementation cost vs Number of Deaths ( W/0 Recovery)
####################

# fig6X, ax6X = plt.subplots(nrows=1, ncols=1 ,sharey = True, sharex = True, figsize=(8, 5), gridspec_kw={'width_ratios': [1]})


# sc = ax6X.scatter( paretoPts_PEPM[:,2], paretoPts_PEPM[:,5], label='PEPM')

# sc = ax6X.scatter( paretoPts_FENM[:,2], paretoPts_FENM[:,5], label='FENM')                

# sc = ax6X.scatter( paretoPts_NEFM[:,2], paretoPts_NEFM[:,5], label='NEFM')

# ax6X.set_xlabel('Total number of deaths')
# ax6X.set_ylabel('Implementation Cost in \$USD')

# fig6X.suptitle('Cost vs Deaths Pareto', wrap = True, fontsize=16)

# ax6X.legend()
# plt.tight_layout()
# plt.subplots_adjust(top = 0.9,wspace=0.02)


#############
# COLLECT INFORMATION
#############

idxDeathSort_PEPM = np.argsort(paretoPts_PEPM[:,2])
idxDeathSort_FENM = np.argsort(paretoPts_FENM[:,2])
idxDeathSort_NEFM = np.argsort(paretoPts_NEFM[:,2])

deathSort_PEPM = paretoPts_PEPM[:,2][idxDeathSort_PEPM]
deathSort_FENM = paretoPts_FENM[:,2][idxDeathSort_FENM]
deathSort_NEFM = paretoPts_NEFM[:,2][idxDeathSort_NEFM]

##############
### Figure 7: Ratio of cost vs total deaths
#############

fig7, (ax71,ax72,ax73) = plt.subplots(nrows=1, ncols=3 ,sharey = True, sharex = True, figsize=(12, 4), gridspec_kw={'width_ratios': [1,1,1]})



enumSort_PEPM = paretoPts_PEPM[:,12][idxDeathSort_PEPM]
mcSort_PEPM = paretoPts_PEPM[:,13][idxDeathSort_PEPM]
localSort_PEPM = paretoPts_PEPM[:,14][idxDeathSort_PEPM]

ax71.plot( deathSort_PEPM , enumSort_PEPM / localSort_PEPM , label ='Enum' , linestyle = '-',linewidth = 2 , color = 'b')
ax71.plot( deathSort_PEPM , mcSort_PEPM / localSort_PEPM , label = 'Mc' , linestyle ='-', linewidth = 2, color = 'g')
ax71.plot( deathSort_PEPM , localSort_PEPM / localSort_PEPM , label ='Local Opt', linestyle =':',linewidth = 2, color = '#FFA500')



enumSort_FENM = paretoPts_FENM[:,12][idxDeathSort_FENM]
localSort_FENM = paretoPts_FENM[:,14][idxDeathSort_FENM]

ax72.plot( deathSort_FENM , enumSort_FENM / localSort_FENM , label ='Enum' , linestyle = '-',linewidth = 2 , color = 'b')
ax72.plot( deathSort_FENM , localSort_FENM / localSort_FENM , label ='Local Opt', linestyle =':',linewidth = 2, color = '#FFA500')


enumSort_NEFM = paretoPts_NEFM[:,12][idxDeathSort_NEFM]
mcSort_NEFM = paretoPts_NEFM[:,13][idxDeathSort_NEFM]
localSort_NEFM = paretoPts_NEFM[:,14][idxDeathSort_NEFM]

ax73.plot( deathSort_NEFM , mcSort_NEFM / localSort_NEFM , label = 'Mc' , linestyle ='-', linewidth = 2, color = 'g')
ax73.plot( deathSort_NEFM , localSort_NEFM / localSort_NEFM , label ='Local Opt', linestyle =':',linewidth = 2, color = '#FFA500')


ax71.set_title('PEPM')
ax72.set_title('FENM')
ax73.set_title('NEFM')

ax71.set_xlabel('Total number of deaths')
ax72.set_xlabel('Total number of deaths')
ax73.set_xlabel('Total number of deaths')
ax71.set_ylabel('Cost ratio')

ax71.legend()
ax72.legend()
ax73.legend()

plt.tight_layout()
plt.subplots_adjust(top = 0.83,wspace=0.02)



#############
# COLLECT INFORMATION
#############


uCtrlMx_PEPM = np.zeros( ( len(paretoPts_PEPM[:,2]) , tf, 4 ) )
uCtrlMx_FENM = np.zeros( ( len(paretoPts_FENM[:,2]) , tf, 4 ) )
uCtrlMx_NEFM = np.zeros( ( len(paretoPts_NEFM[:,2]) , tf, 4 ) )

for ii in range( len(paretoPts_PEPM) ):
    ee,gg = paretoPts_PEPM[ii,10:12]
    uCtrlMx_PEPM[ii,:,:] = plotData_PEPM['uCtrl'][ int(ee), int(gg) , : ,:]

for ii in range( len(paretoPts_FENM) ):
    ee,gg = paretoPts_FENM[ii,10:12]
    uCtrlMx_FENM[ii,:,:] = plotData_FENM['uCtrl'][ int(ee), int(gg) , : ,:]

for ii in range( len(paretoPts_NEFM) ):
    ee,gg = paretoPts_NEFM[ii,10:12]
    uCtrlMx_NEFM[ii,:,:] = plotData_NEFM['uCtrl'][ int(ee), int(gg) , : ,:]


##############
### Figure 8: Testing Control Measures
#############

uBox = [0,tf, np.min(deathSort_PEPM), np.max(deathSort_PEPM)]

fig8 ,  ( (ax81,ax82,ax83),(ax84,ax85,ax86) )   = plt.subplots(nrows=2, ncols=3 ,sharey = True, sharex = True,figsize=(12, 5))    

ax81.imshow(uCtrlMx_PEPM[idxDeathSort_PEPM][:,:,0], 
            origin='lower',
            interpolation='gaussian',
            aspect = 'auto',
            vmin = 0, vmax = uMax[0],
            extent=uBox)

#ax81.set_xlabel('Day')
ax81.set_ylabel('Testing (Low Risk)\n\nTotal number of deaths',wrap = True)
ax81.set_title('PEPM')

ax82.imshow(uCtrlMx_FENM[idxDeathSort_FENM][:,:,0], 
            origin='lower',
            interpolation='gaussian',
            aspect = 'auto',
            vmin = 0, vmax = uMax[0],
            extent=uBox)

#ax82.set_xlabel('Day')
ax82.set_title('FENM')

ax83.imshow(uCtrlMx_NEFM[idxDeathSort_NEFM][:,:,0], 
            origin='lower',
            interpolation='gaussian',
            aspect = 'auto',
            vmin = 0, vmax = uMax[0],
            extent=uBox)

#ax83.set_xlabel('Day')
ax83.set_title('NEFM')

ax84.imshow(uCtrlMx_PEPM[idxDeathSort_PEPM][:,:,1], 
            origin='lower',
            interpolation='gaussian',
            aspect = 'auto',
            vmin = 0, vmax = uMax[1],
            extent=uBox)

ax84.set_xlabel('Day')
ax84.set_ylabel('Testing (High Risk)\n\nTotal number of deaths')
#ax84.set_title('Testing (Low Risk)')

ax85.imshow(uCtrlMx_FENM[idxDeathSort_FENM][:,:,1], 
            origin='lower',
            interpolation='gaussian',
            aspect = 'auto',
            vmin = 0, vmax = uMax[1],
            extent=uBox)

ax85.set_xlabel('Day')
#ax85.set_title('Testing (Low Risk)')

im8 = ax86.imshow(uCtrlMx_NEFM[idxDeathSort_NEFM][:,:,1], 
            origin='lower',
            interpolation='gaussian',
            aspect = 'auto',
            vmin = 0, vmax = uMax[1],
            extent=uBox)

ax86.set_xlabel('Day')
#ax86.set_title('Testing (Low Risk)')


plt.tight_layout()
plt.subplots_adjust(wspace=0.02,right=0.88)

cbar_ax = fig8.add_axes([0.90, 0.15, 0.05, 0.7])
cbar = fig8.colorbar(im8, cax=cbar_ax, orientation='vertical')
cbar.set_label('Control Level')


##############
### Figure 9: Social Distancing Control Measures
#############

uBox = [0,tf, np.min(deathSort_PEPM), np.max(deathSort_PEPM)]

fig9 ,  ( (ax91,ax92,ax93),(ax94,ax95,ax96) )   = plt.subplots(nrows=2, ncols=3 ,sharey = True, sharex = True,figsize=(12, 5))    

ax91.imshow(uCtrlMx_PEPM[idxDeathSort_PEPM][:,:,2], 
            origin='lower',
            interpolation='gaussian',
            aspect = 'auto',
            vmin = 0, vmax = uMax[2],
            extent=uBox)

#ax81.set_xlabel('Day')
ax91.set_ylabel('Distancing (Low Risk)\n\nTotal number of deaths',wrap = True)
ax91.set_title('PEPM')

ax92.imshow(uCtrlMx_FENM[idxDeathSort_FENM][:,:,2], 
            origin='lower',
            interpolation='gaussian',
            aspect = 'auto',
            vmin = 0, vmax = uMax[2],
            extent=uBox)

#ax82.set_xlabel('Day')
ax92.set_title('FENM')

ax93.imshow(uCtrlMx_NEFM[idxDeathSort_NEFM][:,:,2], 
            origin='lower',
            interpolation='gaussian',
            aspect = 'auto',
            vmin = 0, vmax = uMax[2],
            extent=uBox)

#ax83.set_xlabel('Day')
ax93.set_title('NEFM')

ax94.imshow(uCtrlMx_PEPM[idxDeathSort_PEPM][:,:,3], 
            origin='lower',
            interpolation='gaussian',
            aspect = 'auto',
            vmin = 0, vmax = uMax[3],
            extent=uBox)

ax94.set_xlabel('Day')
ax94.set_ylabel('Distancing (High Risk)\n\nTotal number of deaths')
#ax84.set_title('Testing (Low Risk)')

ax95.imshow(uCtrlMx_FENM[idxDeathSort_FENM][:,:,3], 
            origin='lower',
            interpolation='gaussian',
            aspect = 'auto',
            vmin = 0, vmax = uMax[3],
            extent=uBox)

ax95.set_xlabel('Day')
#ax85.set_title('Testing (Low Risk)')

im9 = ax96.imshow(uCtrlMx_NEFM[idxDeathSort_NEFM][:,:,3], 
            origin='lower',
            interpolation='gaussian',
            aspect = 'auto',
            vmin = 0, vmax = uMax[3],
            extent=uBox)

ax96.set_xlabel('Day')
#ax86.set_title('Testing (Low Risk)')


plt.tight_layout()
plt.subplots_adjust(wspace=0.02,right=0.88)

cbar_ax = fig9.add_axes([0.90, 0.15, 0.05, 0.7])
cbar = fig9.colorbar(im9, cax=cbar_ax, orientation='vertical')
cbar.set_label('Control Level')






##############
### Figure 10.1: 
#############
fmt = lambda x,args:'{:.1f}'.format(x*100) 

diff1 = (finalCost_FENM - finalCost_NEFM) / finalCost_NEFM
# neg = diff1 < 0
# diff1 = np.log10( np.maximum( np.abs( diff1 ) , 0.09 ) )
# diff1[neg] = -1 * diff1[neg]

diff2 = (finalCost_PEPM - finalCost_NEFM) / finalCost_NEFM
# neg = diff2 < 0
# diff2 = np.log10( np.maximum( np.abs( diff2 ) , 0.09 ) )
# diff2[neg] = -1 * diff2[neg]

#diff3 = (finalCost_FENM - finalCost_PEPM) / finalCost_PEPM
# neg = diff3 < 0
# diff3 = np.log10( np.maximum( np.abs( diff3 ) , 0.09 ) )
# diff3[neg] = -1 * diff3[neg]

#diffVmin =    #np.min( [np.min(diff1), np.min(diff2) ,np.min(diff3)])
#diffVmax =    #np.max( [np.max(diff1), np.max(diff2) ,np.max(diff3)])



fig10, ax101 = plt.subplots(figsize=(6, 4))


im101 = ax101.imshow( diff1  , cmap='RdBu',origin='lower', interpolation='gaussian', aspect='auto',  extent = box, vmin = -np.max(diff1),vmax=np.max(diff1))

ax101.set_xlabel('$Log_{10}$( Non-Immunity cost )')
ax101.set_ylabel('$Log_{10}$( Death cost )')
ax101.set_title('Percentage cost difference between FENM and NEFM')

cBar1 = fig10.colorbar(im101, cmap='RdBu',orientation='vertical', format= ticker.FuncFormatter(fmt) )
cBar1.set_label('100*(FENM - NEFM)/NEFM')

plt.tight_layout()

##############
### Figure 10.2: 
#############

fig11, ax111 = plt.subplots(figsize=(6, 4))
im111 = ax111.imshow( diff2  ,cmap='RdBu', origin='lower', interpolation='gaussian', aspect='auto',  extent = box,vmin = -np.max(diff2),vmax=np.max(diff2))#, vmin = diffVmin,vmax=diffVmax)

ax111.set_xlabel('$Log_{10}$( Non-Immunity cost )')
ax111.set_ylabel('$Log_{10}$( Death cost )')
ax111.set_title('Percentage cost difference between PEPM and NEFM')

cBar1 = fig11.colorbar(im101, cmap='RdBu',orientation='vertical',format= ticker.FuncFormatter(fmt) )
cBar1.set_label('100*(PEPM - NEFM)/ NEFM')

plt.tight_layout()

##############
### Figure 10.3: 
#############




# fig12, ax121 = plt.subplots(figsize=(6, 4))
# im121 = ax121.imshow( diff3  , cmap='RdBu',origin='lower', interpolation='gaussian', aspect='auto',  extent = box,vmin = -np.max(diff3),vmax=np.max(diff3))#, vmin = diffVmin,vmax=diffVmax)

# ax121.set_xlabel('$Log_{10}$( Non-Immunity cost )')
# ax121.set_ylabel('$Log_{10}$( Death cost )')
# ax121.set_title('FENM - PEPM')

# cBar1 = fig12.colorbar(im121, cmap='RdBu',orientation='vertical', format= ticker.FuncFormatter(fmt) )
# cBar1.set_label('100*(FENM - PEPM)/NEFM')

# plt.tight_layout()





plt.show()










#cBar4 = fig4.colorbar(im43, orientation='vertical',fraction=0.046, pad=0.04 )