#!/usr/bin/env python3
import os
import sys
sys.path.append( '..' )
import MotifTable
import pickle
import distributions

import GenomeBindingTable as gbt
import FragExtract as Frag
import ChipSeq
import PCR

import numpy as np
import pandas as pd
import scipy
import matplotlib.pyplot as plt

cbcolors = {'sky blue': (86/255.0,180/255.0,233/255.0)}

def sampleAccessibility( alpha, numSamples=1000 ):
    """
    Sample values of chromatin accessibility across the human genome. 

    Inputs : 
    alpha --- A scalar value that determines the sensitivity with which
    chromatin accessibility varies with DNAse read counts.

    Returns : 
    A set of "numSamples" values that lie between 0 and 1.  
    """
    df = pd.read_table( 'data/dnase.sample.bed', sep="\t", names=['chr','start','end','read_count'])
    #Divide read counts by 150 to get a DNAse density.
    readCounts = df['read_count'].values/150

    accessibility = np.exp( -alpha*(readCounts ) )

    uniformSet = np.random.random( numSamples )  
    samples = np.zeros_like( uniformSet )

    accessibility.sort()
    quantiles = np.arange( len( accessibility ) )*1.0/(len(accessibility) - 1)

    idx = 0
    for u in uniformSet:
        accIdx = np.argmin( abs(quantiles - u) ) 
        samples[idx] = accessibility[accIdx]
        idx += 1

    return samples

def generateChIPreplicates( spEnergies, numCells=100000, depth=100, ampRatio=1000, pExt=1.0, pcrCycles=15, bgEnergy=1, chemicalPotential=3, kdeReplicates=100 ):
    """
    This function generates multiple replicates of read count ratios from across the genome. 

    """
    numLocations = len( spEnergies )
    simulatedRatios = np.zeros( (kdeReplicates,numLocations) )
    for replicate in range(kdeReplicates):
        if (replicate+1)%5 == 0:
            print( "{}/{}".format( replicate+1, kdeReplicates ) )

        genome = performChipSeq( spEnergies=spEnergies, bgEnergy=bgEnergy, pExt=pExt,  ampRatio=ampRatio, numCells=numCells, pcrCycles=pcrCycles, depth=depth, chemicalPotential=chemicalPotential )
        simulatedRatios[replicate,:] = genome['ratio'].values

    return simulatedRatios

def samplePosteriorSingleTF( spEnergy, simulatedRatios, spEnergies, thinning=10, nMCitr=10000, numTrials=5, prior='uniform',ampRatio=1000,depth=100, numCells=100000,
                            pExt=1.0, pcrCycles=15, eps=0.2, bgEnergy=1, maxReplicates=-1, chemicalPotential=0,makePlot=False,priorParams=[] ):

    #Fraction of iterations that are to be discarded as the burn-in period of
    #the Markov chain.
    burnin = 0.1

    #Since the burn-in period will be discarded, we run the MCMC for 
    #more iterations to ensure that we return "nMCitr" samples to the
    #user.
    totalNumIterations = np.int( (1 + burnin) * nMCitr * thinning )

    #Number of binding locations in the genome.
    numLocations = len( spEnergies )

    kdeList = []
    print("Generating KDEs for different binding probabilities.")
    undefMin = max(spEnergies) + 1
    undefLoc = 0
    obsIdx = np.argmin( abs(spEnergies - spEnergy) )    
    for location in range(numLocations):
        #"mask" contains genomic locations where the read count ratios may be a NaN
        #or an Inf. An Inf occurs when the read count ratio at a location in the input sample 
        #is zero, while a NaN occurs if the read count ratio at a location is zero
        #in both ChIP and input samples. Since these values cannot be used for 
        #kernel density estimation, we keep track of them and discard them. 
        mask = np.isinf( simulatedRatios[:,location] ) | np.isnan( simulatedRatios[:,location] )
        values = simulatedRatios[~mask,location]
        try:
            toAppend = scipy.stats.gaussian_kde( values )
            kdeList.append( toAppend )
        except Exception as info:
            #This exception is triggered when there is only a single
            #read count ratio value present at a location. This is generally because
            #the binding energy of the location is very high. We store the
            #largest energy at which this occurs, and ensure that the MCMC sampler
            #does not generate a candidate value higher than this.
            energy = spEnergies[location]
            undefLoc = location
            if energy < undefMin:
                undefMin = energy

            kdeList.append( None )

    #Since we use a pre-computed read count ratio table, the exact binding
    #energy from which we simulate additional read count ratios may not be
    #present. For instance, there may be no location with a binding energy of
    #exactly 2k_BT in our pre-run simulations. So we pick the genomic location
    #whose pre-assigned binding energy is closest to 2k_BT, say 2.03k_BT, and change it to
    #2k_BT. This does not alter the read count ratios from second genomic locations.
    obsIdx = np.argmin( abs(spEnergies - spEnergy) )    
    spEnergies[obsIdx] = spEnergy

    #We sort binding energies across the genome in order to perform linear
    #interpolation. 
    sortedIdxes = np.argsort( spEnergies )
    sortedEnergies = spEnergies[sortedIdxes]
    sortedKdes = np.array( kdeList )[sortedIdxes]

    trial = 1
    redoPrevious = False
    posteriorIntervalSet = []
    while trial <= numTrials:
        numReplicatesNeeded = 1
        numAvg = 0

        savedSamples = []
        posteriorIntervals = []

        upperLimit = sortedEnergies[-1]
        lowerLimit = sortedEnergies[0]
        while True:
            if not redoPrevious:
                genome = performChipSeq( spEnergies=spEnergies, bgEnergy=bgEnergy, pExt=pExt,  ampRatio=ampRatio, numCells=numCells, pcrCycles=pcrCycles, depth=depth, chemicalPotential=chemicalPotential )
                readRatio = genome.loc[obsIdx,'ratio']
            else:
                redoPrevious = False

            trajectory = np.zeros( totalNumIterations)
            samples = np.zeros( totalNumIterations )
            logAcceptedProb = np.log( np.random.random( totalNumIterations ) )
            spEnergyInit = max(spEnergies) * 0.99 * np.random.random( )
                 
            print("---------------------------------------------------------------")
            print("Observation {} : {}".format( numReplicatesNeeded, readRatio ) )
            print("Starting guess : {}".format( spEnergyInit ) )

            if numReplicatesNeeded > 1:
                logPriorCache = np.copy( logPosteriorCache )
            else:
                #This corresponds to a uniform prior between 0 and 1.
                if prior == 'uniform':
                    logPriorCache = np.ones( numLocations ) * np.log( 1.0/(max(spEnergies)-min(spEnergies)) )
                elif prior == 'powerLaw':
                    #This is a power law distribution as a prior, which is also
                    #the distribution used to sampled binding energies across the genome.
                    k, minX, maxX = priorParams
                    x = np.linspace( priorParams[1]+1e-6, priorParams[2]-1e-6, numLocations )
                    logPriorCache = np.log( distributions.truncPowerLawPdf( x, k, minX, maxX) )

            logLikelihoodCache = []
            bannedList = []
            for idx in range(numLocations):
                if sortedKdes[idx] is None:
                    logLikelihoodCache.append( np.nan )
                    bannedList.append( idx )
                else:
                    logLikelihood = sortedKdes[idx].logpdf( readRatio )[0]
                    logLikelihoodCache.append( logLikelihood )

            logPosteriorCache = []
            for idx in range(numLocations):
                if sortedKdes[idx] is None:
                    logPosteriorCache.append( np.nan )
                else:
                    logPosterior = logLikelihoodCache[idx] + logPriorCache[idx]
                    logPosteriorCache.append( logPosterior )

            idx = np.searchsorted( sortedEnergies, spEnergyInit )
            logProposedPosterior = logPosteriorCache[idx]
            trajectory[0] = logProposedPosterior
            samples[0] = sortedEnergies[idx]
            mcItr = 1

            while mcItr <= totalNumIterations-1:
                proposalEnergy = distributions.truncNorm( a=lowerLimit, b=upperLimit, mu=samples[mcItr-1], sigma=eps, size=1  )[0]

                pdfForward = np.log( distributions.truncNormPdf( proposalEnergy, a=lowerLimit, b=upperLimit, mu=samples[mcItr-1], sigma=eps ) )
                pdfBack = np.log( distributions.truncNormPdf( samples[mcItr-1], a=lowerLimit, b=upperLimit, mu=proposalEnergy, sigma=eps ) )

                idx = np.searchsorted( sortedEnergies, proposalEnergy ) - 1
                slope = (logPosteriorCache[idx+1] - logPosteriorCache[idx])/(sortedEnergies[idx+1] - sortedEnergies[idx])
                interpLogPosterior = logPosteriorCache[idx] + slope * (proposalEnergy - sortedEnergies[idx])
                logProposedPosterior = interpLogPosterior
                if logAcceptedProb[mcItr] < logProposedPosterior + pdfBack - trajectory[mcItr-1] - pdfForward:
                    trajectory[mcItr] = logProposedPosterior
                    samples[mcItr] = proposalEnergy
                else:
                    samples[mcItr] = samples[mcItr-1]
                    trajectory[mcItr] = trajectory[mcItr-1]

                mcItr += 1

            samplesToSave = samples[np.int(burnin * nMCitr * thinning):][::thinning]
            
            if np.unique( samplesToSave ).shape[0] == 1:
                print( "MCMC trapped in poorly sampled range of posterior distribution" )
                print( np.unique( samplesToSave )[0] )
                print( "Restarting MCMC" )
                last = np.unique(samplesToSave)[0]
                if last > spEnergy:     #In this case, lower the upper limit
                    if last < upperLimit:
                        print("Lowering the upper limit to {} from {}".format( last, upperLimit ) )
                        upperLimit = last
                else: #In this case, raise the lower limit
                    if last > lowerLimit:
                        print("Increasing the lower limit to {} from {}".format( last, lowerLimit ) )
                        lowerLimit = last
            
                print("Value of upper limit : {}".format( upperLimit ) )
                redoPrevious = True
                continue
            else:
                logPriorKde = scipy.stats.gaussian_kde( samplesToSave )
                energyEstLower = np.percentile( samplesToSave, 2.5 )
                energyEstHigher = np.percentile( samplesToSave, 97.5 )
                print("95% credible interval : ({},{})".format( energyEstLower, energyEstHigher ))
                print("95% credible interval width : {}".format( energyEstHigher - energyEstLower ))
                print("--------------")

                if makePlot:
                    savedSamples.append( samplesToSave )

                posteriorIntervals.append( [energyEstLower,energyEstHigher] )

                if spEnergy < energyEstLower or spEnergy > energyEstHigher: 
                    trial -= 1
                    print("Repeating trial since 95% credible interval does not include {}".format( spEnergy ) )
                    break
                elif numReplicatesNeeded == maxReplicates:
                    posteriorIntervalSet.append( posteriorIntervals )
                    break
                else:
                    numReplicatesNeeded += 1

        idx = 0
        trial += 1
    
        if not makePlot:
            continue

        #This plots the posterior distribution after four replicates of read count ratios are simulated.
        plt.figure()
        axAllReps = plt.subplot2grid( (numReplicatesNeeded,2), (0,1), rowspan=numReplicatesNeeded )
        axesList = []
        for rep in range(numReplicatesNeeded):
            ax = plt.subplot2grid( (numReplicatesNeeded,2), (rep,0))
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            axesList.append( ax )

        for ax in [axAllReps]:
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)

        for samplesToSave in savedSamples:
            pd.Series( samplesToSave ).plot( kind='kde', ax=axesList[idx], label=str(numReplicatesNeeded), color='black' )
            yLow, yHigh = axesList[idx].get_ylim()
            axesList[idx].plot( posteriorIntervals[idx], [yLow+0.5,yLow+0.5], color='blue', lw=2 )
            axesList[idx].plot( [spEnergy,spEnergy], [yLow,yHigh], color='black', ls='dashed', lw=2 )
            axesList[idx].text( 0.9, 0.9, 'After replicate {}'.format( idx+1 ), ha='right', transform=axesList[idx].transAxes )
            axesList[idx].set_ylim(bottom=0)
            axesList[idx].set_ylabel('')
            if idx < numReplicatesNeeded-1:
                axesList[idx].set_xticklabels([])
            idx += 1
        axesList[idx-1].set_xlabel('Binding energy estimate')

        plt.sca(axAllReps)

    return posteriorIntervalSet

def makeArray( val, N ):
    if type( val ) is list:
        val = np.array( val )
    elif type(val) not in [np.ndarray]:
        val = np.ones( N ) * val

    return val

def performChipSeq( numCells=1000, nReads=100, ampRatio=1000, pExt=1.0, pcrCycles=15,  chromAccessibility=[] ):
    """
    This function combines the GenomeBindingTable,FragExtract,PCR and ChIPseq
    classes together and returns a dataframe that contains read counts at
    each genomic location. Note that this method differs slightly from the
    performChipSeq() function that is part of Animate in that some locations
    can be labelled as false positive binding sites.
    """
    table = gbt.GenomeBindingTable( numCells, chromAccessibility=chromAccessibility  )
    N = len( chromAccessibility )

    pExt = makeArray( pExt, N )
    fragExtract = Frag.FragExtract( pExt, table )
    pAmp = np.power( ampRatio, 1.0/pcrCycles ) - 1 
    pAmp = np.round( makeArray( pAmp, N ), 2 )
    pcrObj = ChipSeq.PCR( pcrCycles, pAmp )

    numReads = makeArray( nReads, numCells )
    chipSeq = ChipSeq.ChipSeq( fragExtract, pcrObj, nReads=numReads )
    
    return chipSeq.uniqueMat

def posteriorEstimate( readMat, priorParams, posteriorDist='beta', posteriorAlpha=0.95 ):
    numCells = readMat.shape[0]

    if posteriorDist in ['beta','gamma-sc','gamma-bulk']:
        alphaPrior, betaPrior = priorParams
        totalReadsPerLocation = readMat.sum(axis=0)

    if posteriorDist == 'beta':
        alphaUpdated = alphaPrior + totalReadsPerLocation
        betaUpdated = betaPrior + numCells - totalReadsPerLocation
        intervals = scipy.stats.beta.interval( posteriorAlpha, alphaUpdated, betaUpdated )
        medians = scipy.stats.beta.median( alphaUpdated, betaUpdated )
    elif posteriorDist == 'gamma-sc':
        alphaUpdated = alphaPrior + totalReadsPerLocation
        betaUpdated = betaPrior/(numCells * betaPrior + 1)        #Reciprocal needed when using scipy.stats.gamma
        intervals = scipy.stats.gamma.interval( posteriorAlpha, alphaUpdated, loc=0, scale=betaUpdated )
        medians = scipy.stats.gamma.median( alphaUpdated, loc=0, scale=betaUpdated )
    elif posteriorDist == 'gamma-bulk':
        alphaUpdated = alphaPrior + totalReadsPerLocation
        betaUpdated = betaPrior/(betaPrior + 1)        #Reciprocal needed when using scipy.stats.gamma
        intervals = scipy.stats.gamma.interval( posteriorAlpha, alphaUpdated, loc=0, scale=betaUpdated )
        medians = scipy.stats.gamma.median( alphaUpdated, loc=0, scale=betaUpdated )

    return [intervals,medians]
     
def main():
    numLocations = 300
    accessibility = sampleAccessibility( 1, numSamples=numLocations )
    accessibility.sort()
    accessibility = accessibility[::-1]
    pExt = 1

    #readsByCell = 50
    numCells = 3000
    readsByCell = 10 + 390 * np.random.random( size=numCells )
    readsByCell = readsByCell.astype(np.uint)

    #uniqueMat is a matrix of dimension (# cells x # locations)
    uniqueMat = performChipSeq( chromAccessibility=accessibility, numCells=numCells, ampRatio=500, pExt=pExt, nReads=readsByCell )

if __name__ == "__main__":
    main()
