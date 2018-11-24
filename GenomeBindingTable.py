from scipy.stats import binom, beta, norm
from scipy.optimize import fsolve
import numpy as np
import pandas as pd

class GenomeBindingTable:
    def __init__( self, numCells, chromAccessibility ):
        """
        The GenomeBindingTable class stores the number of bound fragments based
        on the number of bound fragments in ChIP and input samples at each
        genomic location. 

        To initialize a GenomeBindingTable object, the following non-keyword
        arguments are required inputs.
        1) sequences --- Can be an empty list. Otherwise, a list of N sequences
        where N is the number of genomic locations. The sequences need not be
        equal in length. 
        2) spEnergies --- Binding energies of the target TF at each genomic
        location in the ChIP sample.
        3) bgEnergy --- Binding energy of each location in the input sample.
        4) chemicalPotential --- A value that is proportional the concentration
        of the TF.  
        5) numCells --- Number of cells to be employed in the ChIP sample.

        The keyword arguments are
        6) unboundEnergy --- This is the binding energy of the unbound state (in
        units of kBT) of a genomic location. By default, this is set to 1.59, so
        that the occupancy of the highest affinity site (i.e. site with zero
        energy) is 0.99. 
        6) controlCellRatio --- This is a fraction that determines the number of cells in the
        ChIP sample that will be employed in the control sample. The default
        value is 1.0 i.e. the same number of cells will be employed in both ChIP
        and input samples.
        7) secondTFspEnergies --- Binding energies of the second TF. 
        8) secondTFchemicalPotential --- Chemical potential of the second TF. 
        9) secondTFintEnergies --- An array that specifies the interaction
        energy between both TFs at each genomic location. Positive values
        indicate a competitive interaction, negative values indicate a
        cooperative interaction and zero indicates no interaction. 
        10) indirectLocations --- An array of location numbers that are
        to be simulated as being indirectly bound.  
        11) chromAccessibility -- An array of values that specify the chromatin
        accessibility at each genomic location. The values must lie between
        0 and 1. 

        The function returns a GenomeBindingTable instance containing a
        dataframe called "locations". This dataframe contains N rows and the following
        columns --- 
        1) p_occ_chip --- Probability of occupancy at each location in the ChIP sample
        2) p_occ_bg --- Probability of occucpancy at each location in the
        control sample.
        3) chip_fragments --- Number of fragments bound in the ChIP sample at
        each location 
        4) control_fragments --- Number of fragments bound in the control sample
        at each location.

        Over and above these columns, other columns store the data that was
        input to the GenomeBindingTable object
        5) name --- A number between 1 and N that is assigned to each of the
        N genomic locations.
        6) energy_A --- Binding energy of A that was stored in the spEnergies
        argument.
        7) energy_B --- Binding energy of B that was stored in the
        secondTFspEnergies.
        8) binding --- This can take the following values at each location -
            "direct" (if the first TF binds the location directly or independently of TF B if present), 
            "cooperative" (the first TF binds cooperatively with the second),
            "competitive" (the first TF binds competitively with the second),
            or "indirect" (the first TF binds DNA via the second TF).
        9) int_energy --- The interaction energies stored in the
        secondTFintEnergies argument.
        10) sequence --- The sequence of each genomic location that was stored
        in the "sequences" argument. This is an empty string in each location
        if not sequence was passed.
        """

        #Total number of cells used in the ChIP sample.
        self.numCells = numCells

        self.chromAccessibility = chromAccessibility

        #Total number of binding locations
        self.N = np.int64( len(self.chromAccessibility) )        

        #pTFbound are the probabilities of finding the target TF bound at each
        #genomic location.  pBgBound is the probability of each of these
        #locations being bound in the input sample of the ChIP-seq experiment.
        self.fragmentMat = np.zeros( (self.numCells,self.N), dtype=np.uint8 )
        for row in range(self.fragmentMat.shape[0]):
            self.fragmentMat[row,:] = binom.rvs( 1, self.chromAccessibility )
