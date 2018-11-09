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

def performChipSeq( sequences=[], spEnergies=[], numCells=100000, depth=100, ampRatio=1000, pExt=1.0, pcrCycles=15, bgEnergy=1, chemicalPotential=0, secondTFspEnergies=[], secondTFchemicalPotential=0, chromAccessibility=[], secondTFintEnergies=[], indirectLocations=[], indirectSequences=[], numFPlocations=0 ):
    """
    This function combines the GenomeBindingTable,FragExtract,PCR and ChIPseq
    classes together and returns a dataframe that contains read counts at
    each genomic location. Note that this method differs slightly from the
    performChipSeq() function that is part of Animate in that some locations
    can be labelled as false positive binding sites.
    """
    N = len( spEnergies )
    pExtChip = pExt
    pExtControl = pExt

    numChipReads = N * depth
    numControlReads = N * depth

    bgEnergy = makeArray( bgEnergy, N )
    pAmp = np.round( np.power(ampRatio,1.0/pcrCycles) - 1, 2 )
    pAmp = np.maximum( 0.01, pAmp )
    pAmp = np.minimum( pAmp, 0.99 )

    table = gbt.GenomeBindingTable( sequences, spEnergies,
                                   bgEnergy, chemicalPotential, numCells,
                                   secondTFspEnergies=secondTFspEnergies,
                                   secondTFchemicalPotential=secondTFchemicalPotential,
                                   secondTFintEnergies=secondTFintEnergies,
                                   indirectLocations=indirectLocations,
                                   chromAccessibility=chromAccessibility  )

    pExtControl = makeArray( pExtControl, N )
    pExtChip = makeArray( pExtChip, N )
    fragExtract = Frag.FragExtract( pExtControl, pExtChip, table )
    pAmp = makeArray( pAmp, N )
    pcrObj = ChipSeq.PCR( pcrCycles, pAmp )

    chipSeq = ChipSeq.ChipSeq( table, fragExtract, pcrObj, 
                              nControlReads=numControlReads,
                              nChipReads=numChipReads )
    
    genome = table.locations.merge( chipSeq.readsTable )
    genome = genome.merge( chipSeq.amplifiedTable )
    genome = genome.merge( fragExtract.extractedTable )
    genome.loc[:,'ratio'] = genome.eval('unique_chip_reads/unique_control_reads')

    if numFPlocations > 0:
        genome.loc[:,'ratio'] = genome.eval('unique_chip_reads/unique_control_reads')
        genome.loc[(N-numFPlocations):,'binding'] = 'false-positive'            

    return genome

def getBaselineMotif( chemicalPotential=3, numLocations=1000, tf='Tye7', energiesSampled=[] ):
    """
    This function is used to compute the baseline motif in a ChIP-seq experiment
    where there is no heterogeneity in extraction and amplification efficiency
    across the genome. 

    Returns :  
    There are three items returned in a list in the output. 
    The first are the binding energies employed in the simulation.
    The second is the binding site sequences.
    The third is the weight matrix, or the baseline motif. 
    """
    motif = MotifTable.MotifTable( tf )

    if len( energiesSampled ) == 0:
        energiesSampled = distributions.truncPowerLawRVs( 0.5, 0, 10, size=numLocations )

    spEnergies, sequences = motif.sampleFromDistribution( energiesSampled )
    pcrCycles = 15
    ampRatio = 1000
    genome = performChipSeq( sequences, spEnergies, chemicalPotential=chemicalPotential, ampRatio=ampRatio )

    learnGenome = genome.sort_values(by='ratio',ascending=False).head(np.int(0.1*numLocations))
    pwmSingle, pcmSingle, pfmRef = findPWM( learnGenome['sequence'].values )
    motifSequences = learnGenome['sequence'].values

    return [spEnergies,sequences,motifSequences,pfmRef]

def findPWM( sequences ):
    #Background frequencies of A, T, G, C from the S. cerevisiae genome.
    bgFrequencies = [3.093e-01,1.907e-01,1.907e-01,3.093e-01] 

    numSeq = len( sequences )
    seqLen = len( sequences[0] )

    seqMat = np.zeros( (numSeq,seqLen), dtype=np.object ) 
    for idx in range(numSeq):
        seqMat[idx,:] = list(sequences[idx])

    pwmMat = np.zeros( (4,seqLen), dtype=np.float64 ) 
    pcmMat = np.zeros_like( pwmMat )
    pfmMat = np.zeros_like( pwmMat )
    for col in range(seqLen):
        row = 0
        for letter in ['A','C','G','T']:
            count = np.sum( seqMat[:,col] == letter )*1.0
            pfmMat[row,col] = (count + bgFrequencies[row] )/( bgFrequencies[row] * ( numSeq + 1 ) )
            pcmMat[row,col] = count
            row += 1
        pfmMat[:,col] = pfmMat[:,col]/(pfmMat[:,col].sum())

    for col in range(seqLen):
        row = 0
        for letter in ['A','C','G','T']:
            pwmMat[row,col] = np.log2( pfmMat[row,col]/bgFrequencies[row] )
            row += 1

    return [pwmMat,pcmMat,pfmMat]

def makeDirectSequences( numLocations, tf, distInfo ):
    numLocations = np.int( numLocations )
    motif = MotifTable.MotifTable( tf )

    directDf = pd.DataFrame()

    if distInfo[0] == 'powerLaw':
        distFunc = distributions.truncPowerLawRVs
        parameter, lower, upper = distInfo[1:4]

    lower = max( lower, motif.minEnergy )
    upper = min( upper, motif.maxEnergy )

    #energiesSampled = -distFunc( parameter, -upper, -lower, size=numLocations )
    energiesSampled = distFunc( parameter, lower, upper, size=numLocations )

    spEnergies, sequences = motif.sampleFromDistribution( energiesSampled )
    directDf.loc[:,'sequence'] = sequences
    directDf.loc[:,'sp_energy_1'] = spEnergies
    directDf.loc[:,'mode'] = 'direct'

    return [directDf, motif]

def makeCoopSequences( numCoopLocations, tfPair, distInfoPair ):
    coopDf = pd.DataFrame({})
    numCoopLocations = np.int( numCoopLocations )

    idx = 1
    motifList = []
    for (tf,distInfo) in zip(tfPair,distInfoPair):
        motif = MotifTable.MotifTable( tf )
        motifList.append( motif )

        if distInfo[0] == 'powerLaw':
            distFunc = distributions.truncPowerLawRVs
            parameter, lower, upper = distInfo[1:4]
        
        energiesSampled = distFunc( parameter, lower, upper, size=numCoopLocations)
        spEnergies, sequences = motif.sampleFromDistribution( energiesSampled )
        coopDf.loc[:,'sp_energy_{}'.format( idx )] = spEnergies
        coopDf.loc[:,'sequence_{}'.format( idx )] = sequences
        idx += 1

    coopDf.loc[:,'mode'] = 'cooperative'
    return [coopDf, motifList]

def makeBackgroundSequences( numLocations, tf ):
    numLocations = np.int( numLocations )
    motif = MotifTable.MotifTable( tf )
    numRows, numCols = motif.pwm.shape

    seqLen = np.int( max( numRows, numCols ) )

    backgroundDf = pd.DataFrame()
    backgroundDf.loc[:,'sequence'] = MotifTable.makeDinucBackground( numLocations, seqLen=seqLen )
    backgroundDf.loc[:,'energy_A'] = motif.getInformation( backgroundDf['sequence'].values )

    return backgroundDf

def singleTFmain(makePlot=False,prior='uniform',priorParams=[],maxReplicates=5):
    chemicalPotential = 3

    numLocations = 1000
    pExt = distributions.truncNorm( a=0, b=1, mu=0.5, sigma=0.05, size=numLocations )

    distInfo = ['powerLaw', 0.5, 0, 10]

    if distInfo[0] == 'powerLaw':
        parameter, lower, upper = distInfo[1:4]
        spEnergies = distributions.truncPowerLawRVs( parameter, lower, upper, size=numLocations )
    elif distInfo[0] == 'exp':
        parameter, lower, upper = distInfo[1:4]
        spEnergies = -exponentialLaw( parameter, -upper, -lower, size=numLocations )
    elif distInfo[0] == 'truncNorm':
        distFunc = truncNorm
        parameters, lower, upper = distInfo[1:]
        mu = parameters[0]
        sigma = parameters[1]
        spEnergies = distributions.truncNorm( a=lower, b=upper, mu=mu, sigma=sigma, size=numLocations )

    kdeReplicates = 1000

    numTFs = 1
    prefix = ''
    if isinstance( pExt, float ) or isinstance( pExt, np.float64 ):
        if pExt == 1.0:
            prefix = 'ideal-'
            print("Loading ideal ChIP-seq data")

    fileName = '{}simulatedRatios-N{}-K{}-distInfo-{}-numTFs-{}-mu{}.npz'.format( prefix, numLocations, kdeReplicates, str(distInfo), numTFs, chemicalPotential )
        
    if not os.path.isfile( os.path.join( 'data', fileName ) ): 
        simulatedRatios = generateChIPreplicates( spEnergies, kdeReplicates=kdeReplicates, pExt=pExt, chemicalPotential=chemicalPotential )

        np.savez( os.path.join( 'data', fileName ), simulatedRatios=simulatedRatios, spEnergies=spEnergies, pExt=pExt )
    else:
        print("Loaded file {}".format( fileName ) )
        archive = np.load( os.path.join( 'data', fileName ) )
        simulatedRatios = archive['simulatedRatios']
        spEnergies = archive['spEnergies']
        pExt = archive['pExt']

    nMCitr = 10000
    posteriorData = {}

    prefix = ''
    if isinstance( pExt, float ) or isinstance( pExt, np.float64 ):
        if pExt == 1.0:
            prefix = 'ideal-'

    dictFileName = '{}posteriorData-N{}-K{}-distInfo-{}-numTFs-{}-prior{}-mu{}.pickle'.format( prefix, numLocations, kdeReplicates, str(distInfo), numTFs, prior, chemicalPotential )

    if makePlot:
        posteriorIntervalSet = samplePosteriorSingleTF( 2, simulatedRatios, spEnergies,
                                            numTrials=1, chemicalPotential=chemicalPotential,
                                            prior=prior,pExt=pExt, priorParams=priorParams,
                                            nMCitr=nMCitr,makePlot=True, maxReplicates=maxReplicates )
    else:
        numTrials = 100
        for spEnergy in [2,3,4,5,6]:
            posteriorIntervalSet = samplePosteriorSingleTF( spEnergy, simulatedRatios, spEnergies,
                                            numTrials=numTrials,pExt=pExt, prior=prior, priorParams=priorParams,
                                            chemicalPotential=chemicalPotential, nMCitr=nMCitr, maxReplicates=maxReplicates )

            posteriorData[spEnergy] = posteriorIntervalSet

            dictfile = open( os.path.join( 'data', dictFileName ), 'wb' )
            pickle.dump( posteriorData, dictfile )
            dictfile.close()

def main():
    singleTFmain(prior='powerLaw',maxReplicates=5,priorParams=[0.5,0,10])

if __name__ == "__main__":
    main()
