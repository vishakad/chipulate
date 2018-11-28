from PCR import PCR
import scipy
import numpy as np
import pandas as pd

class ChipSeq:
    def __init__( self, fragExtract, pcr, nReads=-1 ):
        """
        The ChipSeq class ties together all the other classes in the simulation.
        The class contain two key dataframes that contain the number of
        amplified fragments (amplifiedTable) and the number of reads
        (readsTable) for each genomic location.
        """
        amplifiedMat = self.pcrAmplify( fragExtract, pcr )
        self.uniqueMat = self.sampleReads( amplifiedMat, nReads )

    def pcrAmplify( self, fragExtract, pcr ):
        """
        This function calls the PCR sampleFromPCRdist() routine to simulate PCR
        amplification of extracted fragments in both ChIP and control samples.

        Arguments : 
        1) fragExtract --- An instance of the FragExtract class that contains
        the total number of extracted fragments in the ChIP (the
        "ext_chip_fragments" column) and control (the "ext_control_fragments")
        column.
        2) pcr --- An instance of the PCR class that is initialized with the PCR
        efficiencies at each genomic location and a specified number of
        amplification cycles.
        3) fragmentSetStr --- A string that can be set to "chip" or "control"
        depending on whether amplification is being simulated on fragments in
        the ChIP or input sample, respectively.
        """
        amplifiedMat = pcr.sampleFromPCRdist( fragExtract.fragmentMat )
        amplifiedMat = np.ndarray.astype(amplifiedMat, dtype=np.uint32)

        return amplifiedMat

    def sampleReads( self, amplifiedMat, nReads ):
        """
        This function samples nReads from the amplified fragments in the ChIP
        and control samples and returns the total and unique number of reads at
        each location.
        """
        numLocations = amplifiedMat.shape[1]
        uniqueMat = np.zeros_like( amplifiedMat )

        for cell in range(amplifiedMat.shape[0]):
            totalAmplified = amplifiedMat[cell,:].sum()
            amplified = amplifiedMat[cell,:]

            if totalAmplified > nReads[cell]:
                #Sample nReads from across all genomic locations such that each
                #amplified fragment has an equal probability of being chosen.
                readSample = hyperGeomSample( amplified, nReads[cell] )
            else:
                print("WARNING : The number of amplified fragments is less than the total read count. This can happen if the extraction efficiency, number of cells, or the PCR efficiency is too low.  Conversely, the total read count is perhaps too high given the other parameters.")
                readSample = amplified

            uniqueMat[cell,:] = readSample > 0

        return uniqueMat

def hyperGeomSample( binCounts, totalDrawSize ):
    """
    Draw samples from a multi-variate hypergeometric distribution.

    Inputs ---
    1) binCounts --- An array of n values M1,M2,...,Mn where Mi is the number of
    objects in the i-th bin.
    2) totalDrawSize --- The total number of objects to be sampled.

    Returns an array of n value that represent the number of objects sampled
    from each bin.
    """
    numBins = len( binCounts )
    sample = np.zeros( numBins, dtype=np.int64 )
    totalCount = np.sum( binCounts )
    sampleSize = 0

    ctr = 0
    idx = 0

    while idx < numBins:
        success = binCounts[idx]
        failure = totalCount - success
        drawSize = totalDrawSize - sampleSize

        if drawSize == 0:
            return sample

        draw = np.random.hypergeometric( success, failure, drawSize )
        if draw > totalDrawSize:
            sample[ctr] = n
            return sample
        else:
            totalCount -= binCounts[idx]
            sample[ctr] = draw
            sampleSize += draw

        ctr += 1

        idx += 1

    return sample

def hyperGeomSample2( binCounts, totalDrawSize ):
    """
    Draw samples from a multi-variate hypergeometric distribution.

    Inputs ---
    1) binCounts --- An array of n values M1,M2,...,Mn where Mi is the number of
    objects in the i-th bin.
    2) totalDrawSize --- The total number of objects to be sampled.

    Returns an array of n value that represent the number of objects sampled
    from each bin.
    """
    totalCount = np.sum( binCounts )
    choices = np.random.choice( totalCount, replace=False, size=totalDrawSize )
    sample = np.zeros_like( binCounts )

    cumBinCounts = np.append( [0], np.cumsum( binCounts ) )
    for idx in range(len(binCounts)):
        lb = cumBinCounts[idx]
        ub = cumBinCounts[idx+1]
        sample[idx] = np.count_nonzero( (choices >= lb) & (choices < ub) )

    return sample
