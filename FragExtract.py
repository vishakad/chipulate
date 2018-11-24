from scipy.stats import binom
from scipy.stats import beta
import numpy as np
import pandas as pd
 
class FragExtract:
    def __init__( self, pExt, boundTable ):
        """
        The FragExtract class contains routines needed for simulating fragment
        extraction from a pool of bound fragments. 
        
        The arguments to the class constructor are three non-keyword quantities ---
        pExtControl, pExtChip and boundTable
        1) pExtControl --- An array of values of size N, where N is the number
        of binding locations, which specifies the extraction efficiency at
        each location for the input/control sample. Each value must lie between 0 and 1. 
        2) pExtChip --- An array of values of the same dimension as pExtControl, 
        which specifies the extraction efficiency at each genomic location in
        the ChIP sample. In general, we could set pExtChip and pExtControl to be different
        values but we keep them equal in all simulations.
        3) boundTable --- An instance of the GenomeBindingTable class. The
        instance should have the control_fragments and chip_fragments variables 
        set to non-zero values. 
        """
        #pExtChip - Antibody efficiency/ChIP extraction efficieny -- What is the probability of a fragment that is bound by the TF
        #of interest being picked up by the antibody in the ChIP sample?
        #Default : 1 
        self.pExt = pExt

        self.fragmentMat = np.zeros_like( boundTable.fragmentMat )
        for row in range(boundTable.fragmentMat.shape[0]):
            self.fragmentMat[row,:] = binom.rvs( boundTable.fragmentMat[row,:], self.pExt )
 
