import numpy as np
from os.path import isfile, join
from scipy.sparse import dok_matrix
from scipy.stats import binom

def makeArray( values, N ):
    """
        This function converts scalar values and Python lists into Numpy arrays
        of size 1 x N. 
    """
    if type( values ) is list:
        values = np.array( values )
    elif type(values) not in [np.ndarray]:
        values = np.ones( N ) * values

    return values

class PCR:
    def __init__( self, nCycles, pAmp, nQuantiles=10000 ):
        """
        The PCR class contains routines for simulating PCR amplification of
        extract fragments from across the genome.  

        The inputs to the PCR constructor are ----
        1) nCycles --- The number of cycles of PCR amplification.
        2) pAmp --- A scalar or an array of N PCR efficiencies, where N is the
        number of locations and each efficiency lies between 0 and 1.
        3) nQuantiles --- This sets the accuracy of the PCR simulation process.
        A larger number provides a better accuracy. The default value is set at
        10000, and setting this parameter below 1000 is not recommended. 

        Increasing pAmp or nCycles slows down the simulation process. 
        """
        self.nCycles = nCycles
        self.nQuantiles = np.int( nQuantiles )
        self.pcrDistCdf = {}
        self.pAmp = pAmp
        self.preComputePCRquantiles( )

    def preComputePCRquantiles( self ):
        """
        This function computes the quantile values of the probability mass
        function of the number of amplified fragments. 

        These quantiles are key to the inverse sampling method that is employed
        in simulating PCR amplification.
        """
        quantileSet = np.linspace( 0, 1, self.nQuantiles ) 
        self.pcrDistQuantiles = {}

        if type(self.pAmp) in [float,np.float64]:
            pAmp = self.pAmp
            #np.cumsum computes the cumutuate density function (CDF) of the
            #amplification distribution. The n-th location in
            #self.pcrDistCdf[str(pAmp] is the probability of observing
            #at most n amplified fragments, starting from a single fragment
            #before PCR.
            self.pcrDistCdf[str(pAmp)] = np.cumsum( self.computePCRdist( pAmp ) )
            #Once the CDF is computed, the quantileSet variable contains
            #possible values of the CDF between 0 and 1 in steps of 1/self.nQuantiles.
            #Basically, q = P[ X < n ] where q is the quantile and P[X <= n] is
            #the CDF of X, where X is the random variable denoting the number of
            #amplified fragments.
            #For each value of q between 0 and 1, the np.searchsorted function
            #finds a value of n such that the above equation is satisfied.
            self.pcrDistQuantiles[pAmp] = np.searchsorted( self.pcrDistCdf[str(pAmp)], quantileSet )+1
        else:
            #If there is heterogeneity in PCR efficiency across the genome, then
            #the PCR amplification distribution needs to be computed for each
            #PCR efficiency in the genome. 
            pAmps = np.unique( self.pAmp )
            for pAmp in pAmps:
                self.pcrDistCdf[str(pAmp)] = np.cumsum( self.computePCRdist( pAmp ) )
                self.pcrDistQuantiles[pAmp] = np.searchsorted( self.pcrDistCdf[str(pAmp)], quantileSet )+1

    def sampleFromPCRdist( self, fragmentCounts, singlePamp=None, returnPerFragment=False ):
        N = len( fragmentCounts )

        #In case there is no PCR heterogeneity in the genome, make a numpy array
        #where each location's PCR efficiency is set to the same value. 
        if singlePamp is None:
            pAmps = makeArray( self.pAmp, N )
        else:
            pAmps = makeArray( singlePamp, N )

        fragmentCounts = makeArray( fragmentCounts, N )
        fragmentCounts = fragmentCounts.astype( np.int )
        amplifiedFragments = np.zeros_like( fragmentCounts )

        totalFragments = fragmentCounts.sum()

        #If returnPerFragment is set to True, then the number of
        #amplified fragments per extracted fragment is returned. 
        #For example, if there are 3 genomic locations with
        #the extracted fragment array being [2,3,4], then there
        #are a total of 2 + 3 + 4 = 9 samples drawn from the PCR distribution.
        #Suppose the PCR efficiencies in these three locations are [0.5,0.8,0.9],
        #and the numbers of amplified fragments obtained from each of the 9
        #extracted fragments are [1243,1522,644,3553,3421,2034,6735,7852,9425].
        #When returnPerFragment is True, the perFragmentAmplified
        #list returned is [ [1243,1522], [644,3553,3421], [2034,6735,7852,9425] ].
        if returnPerFragment:
            perFragmentAmplified = np.zeros( N, dtype=np.object )
        else:
            perFragmentAmplified = []

        idx = 0 
        #As part of inverse sampling, we first draw quantiles between 1 and
        #self.nQuantiles with equal probability of choosing each value. 
        quantilesAll = np.random.randint( self.nQuantiles, size=totalFragments )
        quantileIdx = 0
        for fragmentCount in fragmentCounts:
            #fragmentCount refers to the number of pre-PCR fragments at a single
            #genomic location where the PCR efficiency is pAmps[idx].
            quantiles = quantilesAll[quantileIdx:(quantileIdx+fragmentCount)] 
            pcrCopies = 0

            #We choose the quantiles of the PCR amplification distribution
            #whose PCR efficiency correponds to pAmps[idx]. 
            pcrQuantiles = self.pcrDistQuantiles[pAmps[idx]]

            #Since pcrQuantiles stores the quantiles of the PCR distribution,
            #the "quantiles" array can be used to directly simulate PCR. 
            temp = pcrQuantiles[quantiles]
            if returnPerFragment:
                perFragmentAmplified[idx] = temp 
            amplifiedFragments[idx] = temp.sum()
            idx += 1
            quantileIdx += fragmentCount

        if returnPerFragment:
            return [amplifiedFragments, perFragmentAmplified]
        else:
            return amplifiedFragments

    def computePCRdist( self, pAmp ):
        """
        This function computes the probability mass function of the number of
        amplified fragments present when a single fragment of DNA is amplified
        after self.nCycles of amplification are carried out with a PCR
        efficiency of pAmp.
        
        Inputs :
        pAmp - A scalar value of PCR efficiency that lies between 0 and 1. 

        Returns : 
        A 1 x K array of values, where the n-th position is the probability of
        obtaining n fragments. K, which is 2 ^ (number of cycles), is the
        maximum number of fragments obtained if every fragment gets amplified
        after each cycle of amplification. 

        This distribution obtained for the input PCR efficiency is also written
        to disk, which can be retrieved for later use when amplification is
        simulated. 

        See Methods section in the accompanying manuscript for details on the
        computation method. This method is an implementation of the PCR
        simulation scheme proposed in Kebschull, Justus M., and Anthony M.
        Zador. "Sources of PCR-induced distortions in high-throughput sequencing
        data sets." Nucleic acids research 43.21 (2015): e143-e143.
        """

        possFilePath = join( 'output', 'pcr-data', '{}-{}.npy'.format( self.nCycles, pAmp ) )

        #In case the distribution has been already computed before, retrieve it
        #and return it instead of re-computing it afresh.
        if isfile( possFilePath ):
            copyNumberDist = np.load( possFilePath )
            return copyNumberDist

        #If the distribution has not been computed before, then go ahead and
        #compute it.
        maxCopyNumber = (2 ** self.nCycles) + 1
        copyNumberDist = np.zeros( maxCopyNumber ).reshape( (maxCopyNumber,1) )
        copyNumberDist[0] = 1

        pcrMatrix = dok_matrix( (maxCopyNumber, maxCopyNumber), dtype=np.float32  ) 
        copyNumberValues = range( 1,maxCopyNumber+1 )
        print("Computing PCR distribution for {},{}".format( self.nCycles, pAmp ))
        for column in copyNumberValues:
            pcrMatrix[column-1:,column-1] = binom.pmf( range( maxCopyNumber-column+1 ), column, pAmp ).reshape( (maxCopyNumber-column + 1, 1 ) )

        for idx in range( self.nCycles ):
            copyNumberDist = pcrMatrix * copyNumberDist

        copyNumberDist = copyNumberDist/np.sum( copyNumberDist )
        np.save( possFilePath, copyNumberDist )

        return copyNumberDist
