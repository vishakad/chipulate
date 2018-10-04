import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib import cm
from matplotlib import gridspec
from matplotlib import colors as mplcolors

import CBcm
import numpy as np
import util

rcParams['svg.fonttype'] = 'none'
def makeColours( vals, cmapStr, returnDict=False ):
    colours = np.zeros( (len(vals),3) )
    origVals = np.copy( vals )
    if vals.max() > 1:
        vals = vals/vals.max()
        norm = mplcolors.Normalize( vmin=vals.min(), vmax=vals.max()*0.95 )
    else:
        norm = mplcolors.Normalize( vmin=vals.min(), vmax=vals.max()*0.95 )
    colours = [cm.ScalarMappable( norm=norm, cmap=cmapStr).to_rgba( val ) for val in vals]
    if returnDict:
        uniqueOrigVals = np.unique( origVals )
        uniqueVals = np.unique( vals )
        uniqueColours = [cm.ScalarMappable( norm=norm, cmap=cmapStr).to_rgba( val ) for val in uniqueVals]
        valColourDict = {}

        for (val,colour) in zip(uniqueOrigVals,uniqueColours):
            valColourDict[colour] = val 
        
        return [colours, valColourDict]
    else:
        return colours

def manyErrorBarsPlot( x, yMeans, yStds, lineLevels, cmap=CBcm.CB2cm['YeBl'], xIsMean=False, cbarText='', xAxisText='', yAxisText='', titleText='', ylim=[], fontsize=10, plotCbar=True, ls='solid' ):
    if plotCbar:
        plt.figure()
        gs = mpl.gridspec.GridSpec(1, 2, width_ratios=[8, 1]) 
        axes = [plt.subplot( gsAx ) for gsAx in gs]
    else:
        axes = [plt.gca()]
    
    if isinstance( lineLevels, list):
        lineLevels = np.array( lineLevels )

    colours = makeColours( lineLevels, cmap )
    if xIsMean:
        for idx in range(yMeans.shape[1]):
            axes[0].errorbar( x, yMeans[:,idx], yerr=yStds[:,idx], capthick=2, color=colours[idx], ms=2, marker='.', ls=ls )
    else:
        for idx in range(yMeans.shape[0]):
            axes[0].errorbar( x, yMeans[idx,:], yerr=yStds[idx,:], capthick=2, color=colours[idx], ms=2, marker='.', ls=ls )
        
    axes[0].set_xlabel( xAxisText, fontsize=fontsize )
    axes[0].set_ylabel( yAxisText, fontsize=fontsize )
    axes[0].set_title( titleText )
    if len(ylim) > 0:
        axes[0].set_ylim(ylim)

    if plotCbar:
        #norm = mpl.colors.Normalize( vmin=1, vmax=4 )
        norm = mpl.colors.Normalize( vmin=min(lineLevels), vmax=max(lineLevels) )    
        #norm = mpl.colors.LogNorm( vmin=), vmax=max(lineLevels) )    
    
        cbl = mpl.colorbar.ColorbarBase( axes[1], norm=norm, cmap=cmap, 
                                        orientation='vertical', values=lineLevels )

        tickBoundaries = np.arange( len(lineLevels) )*1.0/len(lineLevels)
        tickBoundaries = np.append( tickBoundaries, 1 )
        tickLocs = []
        for idx in range(len(tickBoundaries)-1):
            tickLocs.append( 0.5 * (tickBoundaries[idx] + tickBoundaries[idx+1]))

        cbl.ax.yaxis.set_ticks( tickLocs )
        cbl.ax.yaxis.set_ticklabels( lineLevels )
        cbl.ax.yaxis.set_label_text( cbarText )
        plt.tight_layout() 
    
    if plotCbar:
        return [axes[0],cbl]
    else:
        return [axes[0],colours]

def getTpr( genome, fprTarget ):
    if fprTarget > 1-1e-6:
        return [1,1,-np.inf]

    readCounts = genome['ratio'].values
    lb = min( readCounts )
    ub = max( readCounts )

    numFP = genome.query('binding == "false-positive"').shape[0]
    numTP = genome.shape[0] - numFP
    #print("Starting search. Target : {}".format( fprTarget ))
    fpr = 1.0
    nItr = 100
    itr = 0
    while abs( fpr - fprTarget ) > 0.01 * fprTarget and itr < nItr:

        threshold = (lb + ub)/2.0
        genome.loc[:,'inference'] = 'no-binding'
        genome.loc[readCounts > threshold,'inference'] = 'binding'

        if threshold < 1e-8:
            fpr = fprTarget
            print( [fpr,fprTarget,threshold])
            break

        fpr = genome.query( '(inference == "binding") & (binding == "false-positive")').shape[0]*1.0/numFP
        if fpr > fprTarget:
            lb = threshold
        else:
            ub = threshold

        itr += 1

    tpr = genome.query( '(inference == "binding") & (binding != "false-positive")').shape[0]*1.0/numTP
    return [fpr,tpr,threshold]

def makeROC( genome, color, ax=None, ls='solid', fprTargets=[], label="" ):
    if len( fprTargets ) == 0:
        fprTargets = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]

    fprVals = []
    tprVals = []
    
    readCountThresholds = []
    for fprTarget in fprTargets:
        if abs( fprTarget ) < 1e-6:
            fpr = 0
            tpr = 0
            readCountThreshold = max( genome['ratio'].values )
        else:
            fpr, tpr, readCountThreshold = getTpr( genome, fprTarget )

        fprVals.append( fpr ) 
        tprVals.append( tpr )
        readCountThresholds.append( readCountThreshold )

    itr = 0
    if ax is None:
        fig = plt.figure()
        ax = plt.gca()

    auROC = util.findAuROC( genome )

    if label == "":
        ax.plot( fprVals, tprVals, color=color, lw=0.75, ls=ls, label='{:.2f}'.format(auROC), marker='o' , markersize=2 )
    else:
        ax.plot( fprVals, tprVals, color=color, lw=0.75, ls=ls, label=label + ' ({:.2f})'.format(auROC), marker='o', markersize=2 )

    ax.set_xlim(0,1)
    ax.grid('on',ls='dashed')
    ax.legend(fontsize=8)

    return [ax,auROC,readCountThresholds]
