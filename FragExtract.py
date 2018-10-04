from scipy.stats import binom
from scipy.stats import beta
import numpy as np
import pandas as pd
 
class FragExtract:
    def __init__( self, pExtControl, pExtChip, boundTable ):
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
        self.pExtChip = pExtChip
 
        #pExtControl - Probability of extracting a fragment from the control sample. 
        self.pExtControl = pExtControl
 
        #extractedTable is a pandas dataframe that stores the number of 
        #extracted fragments in the ChIP and input samples at each genomic location.
        self.extractedTable = pd.DataFrame()
        self.extractedTable.loc[:,'name'] = boundTable.locations['name'].values
             
        self.N = boundTable.locations.shape[0]
 
        #The extracted fragments are binomially sampled from bound fragments.
        extControlFragments = binom.rvs( boundTable.locations['control_fragments'], self.pExtControl, size=self.N )

        self.extractedTable.loc[:,'ext_control_fragments'] = extControlFragments
 
        extChipFragments = binom.rvs( boundTable.locations['chip_fragments'], self.pExtChip, size=self.N )
            
        self.extractedTable.loc[:,'ext_chip_fragments'] = extChipFragments
