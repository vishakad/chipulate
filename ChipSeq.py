from PCR import PCR
import scipy
import numpy as np
import pandas as pd

class ChipSeq:
    def __init__( self, genomeBindingTable, fragExtract, pcr, nChipReads=-1, nControlReads=-1, generateIntervals=True ):
        """
        The ChipSeq class ties together all the other classes in the simulation.
        The class contain two key dataframes that contain the number of
        amplified fragments (amplifiedTable) and the number of reads
        (readsTable) for each genomic location.
        """
        self.amplifiedTable = pd.DataFrame( {'name' : fragExtract.extractedTable['name'].values} ) 
        self.readsTable = pd.DataFrame( {'name' : fragExtract.extractedTable['name'].values} ) 
        self.controlFragmentLocationMatrix = []
        self.chipFragmentLocationMatrix = []
        self.controlReadsLocationMatrix = []
        self.chipReadsLocationMatrix = []
        self.perControlAmplified = []
        self.perChipAmplified = []
        self.bindingTable = genomeBindingTable

        fragExtract = self.downsampleControl(fragExtract)

        self.amplifiedTable.loc[:,'amp_control_fragments'], self.perControlAmplified = self.pcrAmplify( fragExtract, pcr, 'control' )
        self.readsTable.loc[:,'unique_control_reads'], self.readsTable.loc[:,'control_reads'], self.controlFragmentNumbers = self.sampleReads( fragExtract, 'control', nControlReads, generateIntervals )

        self.amplifiedTable.loc[:,'amp_chip_fragments'], self.perChipAmplified = self.pcrAmplify( fragExtract, pcr, 'chip' )
        self.readsTable.loc[:,'unique_chip_reads'], self.readsTable.loc[:,'chip_reads'], self.chipFragmentNumbers = self.sampleReads( fragExtract, 'chip', nChipReads, generateIntervals )

    def downsampleControl( self, fragExtract ):
        """
        This function downsamples either the number of control fragments or number of
        ChIP fragments to ensure that their total numbers are equal.
        """
        #If there are more fragments in the input sample than in the ChIP sample
        if fragExtract.extractedTable['ext_control_fragments'].sum() > fragExtract.extractedTable['ext_chip_fragments'].sum():
            downsampling = fragExtract.extractedTable['ext_chip_fragments'].sum()*1.0/fragExtract.extractedTable['ext_control_fragments'].sum()
            extControlFragments = scipy.stats.binom.rvs( fragExtract.extractedTable['ext_control_fragments'], downsampling, size=fragExtract.extractedTable.shape[0] )

            fragExtract.extractedTable.loc[:,'ext_control_fragments'] = extControlFragments
        else:
            #If there are more fragments in the ChIP sample than the input sample.
            downsampling = fragExtract.extractedTable['ext_control_fragments'].sum()*1.0/fragExtract.extractedTable['ext_chip_fragments'].sum()
            extChipFragments = scipy.stats.binom.rvs( fragExtract.extractedTable['ext_chip_fragments'], downsampling, size=fragExtract.extractedTable.shape[0] )

            if len(extChipFragments) < fragExtract.extractedTable.shape[0]:
                extChipFragments = np.append( extChipFragments, np.zeros( fragExtract.extractedTable.shape[0] - len(extChipFragments) ) )

            fragExtract.extractedTable.loc[:,'ext_chip_fragments'] = extChipFragments

        return fragExtract

    def pcrAmplify( self, fragExtract, pcr, fragmentSetStr ):
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
        fragmentCounts = fragExtract.extractedTable['ext_{}_fragments'.format( fragmentSetStr )]
        numLocations = fragExtract.extractedTable.shape[0]

        amplified, perFragmentAmplified = pcr.sampleFromPCRdist( fragmentCounts, returnPerFragment=True )
        amplified = np.ndarray.astype(amplified, dtype=np.int64)

        return [amplified, perFragmentAmplified]

    def sampleReads( self, fragExtract, fragmentSetStr, nReads, generateIntervals ):
        """
        This function samples nReads from the amplified fragments in the ChIP
        and control samples and returns the total and unique number of reads at
        each location.
        """
        if fragmentSetStr == 'control':
            perFragmentAmplified = self.perControlAmplified
        else:
            perFragmentAmplified = self.perChipAmplified

        amplified = self.amplifiedTable['amp_{}_fragments'.format( fragmentSetStr )].values
        amplified = np.ndarray.astype(amplified, dtype=np.int64)
        unamplified = fragExtract.extractedTable['ext_{}_fragments'.format( fragmentSetStr )]

        totalAmplified = amplified.sum()
        N = len( amplified )

        if totalAmplified > nReads:
            #Sample nReads from across all genomic locations such that each
            #amplified fragment has an equal probability of being chosen.
            readSample = hyperGeomSample( amplified, nReads )
        else:
            print("WARNING : The number of amplified fragments is less than the total read count. This can happen if the extraction efficiency, number of cells, or the PCR efficiency is too low.  Conversely, the total read count is perhaps too high given the other parameters.")
            readSample = amplified

        uniques = np.zeros( N, dtype=np.int )
        duplicates = np.zeros( N, dtype=np.int )

        if generateIntervals:
            readsToDuplicate = np.zeros( int(nReads), dtype=np.int )
        else:
            readsToDuplicate = []

        #See the Methods section in the manuscript for details on how reads 
        #are sampled from the pool of amplified fragments.
        uniqueReadIdx = 0
        for i in range( N ):
            if readSample[i] > 0:
                locsToChoose = hyperGeomSample( perFragmentAmplified[i], readSample[i] )
                mask = locsToChoose >= 1
                uniques[i] = np.sum( mask, dtype=np.int )
                if generateIntervals:
                    readsToDuplicate[uniqueReadIdx:(uniqueReadIdx+uniques[i])] = locsToChoose[mask]
            else:
                uniques[i] = 0

            uniqueReadIdx += uniques[i]

        if generateIntervals:
            readsToDuplicate = readsToDuplicate[:uniqueReadIdx]

        uniques = np.ndarray.astype(uniques,dtype=np.int64)

        if fragmentSetStr == 'chip':
            num = self.amplifiedTable.shape[0]
            uniques = np.append( uniques, np.zeros( num - len(uniques), dtype=np.int64 ) )
            readSample = np.append( readSample, np.zeros( num - len(readSample), dtype=np.int64 ) )

        return [uniques,readSample,readsToDuplicate] 

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
