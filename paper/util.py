import numpy as np
import sklearn.metrics

def computeFidelity( genome, pSpColumn='energy_A', ratioColumn='ratio', minPairs=1000, quantiles=[] ):
    """
    This function computes the fidelity of a ChIP-seq experiment.  
    """
    qIdx = 0
    if len( quantiles ) == 0:
        quantiles = [[0,0.25],[0.25,0.5],[0.5,0.75],[0.75,1.0],[0,1.0]]

    fidelityValues = np.zeros( len(quantiles) )

    if ratioColumn == 'ratio':
        maxRatio = genome.loc[~np.isinf(genome[ratioColumn]) & ~np.isnan( genome[ratioColumn] ),ratioColumn].max()
        locs = genome['unique_control_reads'] == 0
        genome.loc[locs,ratioColumn] = maxRatio + 0.01
        locs = genome['unique_chip_reads'] == 0
        genome.loc[locs,ratioColumn] = 0

    for q in quantiles:
        qLow, qHigh = q
        qLow, qHigh = genome[ratioColumn].quantile( [qLow,qHigh] )
        genomeCopy = genome.query('{} >= {} & {} <= {}'.format( ratioColumn, qLow, ratioColumn, qHigh ) ).copy()
            
        ratioCounts = genomeCopy[ratioColumn].values
        pSpValues = genomeCopy[pSpColumn].values

        numPairsSampled = minPairs * 2

        sampledPairs = np.array([])
        while True:
            pairs = np.random.randint( len(ratioCounts), size=(np.int64(numPairsSampled),2) )
            pairsToRedo = (pairs[:,0] == pairs[:,1])

            if np.sum( pairsToRedo ) > 0:
                signs = ( np.random.random( size=np.sum(pairsToRedo)) > 0.5 ).astype( np.int )
                pairs[pairsToRedo,1] = (pairs[pairsToRedo,0] + np.power(-1,signs) * np.random.randint( low=1, high=100, size=np.sum( pairsToRedo ))) % len(ratioCounts)

            toShuffle = ratioCounts[pairs[:,1]] <= ratioCounts[pairs[:,0]]
            pairs[toShuffle,0], pairs[toShuffle,1] = pairs[toShuffle,1], pairs[toShuffle,0]

            left = pairs[:,0]                
            right = pairs[:,1]

            ratiosLeft = ratioCounts[left]
            ratiosRight = ratioCounts[right]

            criterion = (ratiosRight >= 0.9*ratiosLeft)
            if sampledPairs.shape[0] == 0:
                sampledPairs = pairs[criterion,:]
            else:
                sampledPairs = np.append( sampledPairs, pairs[criterion,:],axis=0)

            if sampledPairs.shape[0] > minPairs:
                sampledPairs = sampledPairs[:minPairs]
                break
            else:
                numPairsSampled *= 1.5

        left = sampledPairs[:,0]
        right = sampledPairs[:,1]
        pLeft = pSpValues[left]
        pRight = pSpValues[right]
        ratiosLeft = ratioCounts[left]
        ratiosRight = ratioCounts[right]

        #mask = (pLeft >= pRight) & (ratiosLeft <= ratiosRight)
        fidelityMask = np.zeros( minPairs )
        idxes = np.arange( minPairs )
        mask = (pLeft >= pRight) & (ratiosLeft <= ratiosRight)
        fidelityMask[mask] = 1
        undecided = idxes[np.isclose(ratiosLeft, ratiosRight)]
        randomChoice = np.random.random( len(undecided) )
        fidelityMask[undecided] = randomChoice >= 0.5
        fidelityMask[undecided] = randomChoice <= 0.5

        fidelityValues[qIdx] = np.sum( fidelityMask )*1.0/np.sum( ratiosLeft <= ratiosRight )

        qIdx += 1

    return fidelityValues

def findAuROC( genome ):
    """
    This function computes the area under the ROC curve of a ChIP-seq experiment
    where read count ratios are employed to distinguish between bound and
    unbound regions.
    """
    labels = np.zeros(genome.shape[0],dtype=np.bool)
    direct = genome['binding'] == "direct"
    labels[direct] = True
    labels[~direct] = False
    auROC = sklearn.metrics.roc_auc_score( labels, genome['ratio'].values )

    return auROC
