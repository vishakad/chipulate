from scipy.stats import binom, beta, norm
from scipy.optimize import fsolve
import numpy as np
import pandas as pd

class GenomeBindingTable:
    def __init__( self, sequences, spEnergies, bgEnergy, chemicalPotential, numCells, names=[], unboundEnergy=1.59, controlCellRatio=0.1, secondTFspEnergies=[], secondTFchemicalPotential=0, secondTFintEnergies=[], indirectLocations=[], chromAccessibility=[]):
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
        6) names --- Names for each genomic region.
        7) unboundEnergy --- This is the binding energy of the unbound state (in
        units of kBT) of a genomic location. By default, this is set to 1.59, so
        that the occupancy of the highest affinity site (i.e. site with zero
        energy) is 0.99. 
        8) controlCellRatio --- This is a fraction that determines the number of cells in the
        ChIP sample that will be employed in the control sample. The default
        value is 1.0 i.e. the same number of cells will be employed in both ChIP
        and input samples.
        9) secondTFspEnergies --- Binding energies of the second TF. 
        10) secondTFchemicalPotential --- Chemical potential of the second TF. 
        11) secondTFintEnergies --- An array that specifies the interaction
        energy between both TFs at each genomic location. Positive values
        indicate a competitive interaction, negative values indicate a
        cooperative interaction and zero indicates no interaction. 
        12) indirectLocations --- An array of location numbers that are
        to be simulated as being indirectly bound.  
        13) chromAccessibility -- An array of values that specify the chromatin
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
        self.numCells = np.float64( numCells )
        self.chipCells = np.int64( numCells )

        #Binding mismatch energies for the target TF
        self.spEnergies = spEnergies                

        #Total number of binding locations
        self.N = np.int64( len(spEnergies) )        

        #Background binding energy for the control sample
        self.bgEnergy = bgEnergy                    

        #Binding energy E_0 corresponding to the unbound state. By default, this
        #is set to a value where the occupancy probability at the highest
        #affinity site is 0.99
        self.unboundEnergy = unboundEnergy

        #Chemical potential of the target TF. Default value : 0
        self.chemicalPotential = chemicalPotential  

        #Fraction of cells that is used for the control experiment. Default : 0.9
        self.controlCells = np.int64( controlCellRatio * self.numCells ) 

        #Binding energies of the second TF that may be passed in 
        #to the simulation.
        self.secondTFspEnergies = secondTFspEnergies

        #Chemical potential of the second TF
        self.secondTFchemicalPotential = secondTFchemicalPotential

        #Interaction energy of the target TF (A) with the second TF (B)
        self.secondTFintEnergies = secondTFintEnergies

        #Indices of locations that are indirectly bound.
        self.indirectLocations = indirectLocations

        #Each location is assigned a name, which is just a number between 1 and N,
        #where N is the number of binding locations. This is the column used
        #to join entries with second tables from the fragment extraction, 
        #PCR amplification and sequencing processes.
        self.locations = pd.DataFrame( columns=['name'] )
        if len( names ) == 0:
            self.locations.loc[:,'name'] = ['region_' + str(idx) for idx in range( 1, self.N+1 )]
        else:
            self.locations.loc[:,'name'] = names

        #Binding energies of the TF A at each location.
        self.locations.loc[:,'energy_A'] = spEnergies

        #By default, all bound locations are assumed to be directly bound, and
        #are assigned the label 'direct'.
        self.locations.loc[:,'binding'] = 'direct'

        if len(sequences) > 0:
            self.locations.loc[:,'sequence'] = sequences
        else:
            self.locations.loc[:,'sequence'] = [""]*self.N

        #Chromatin accessibility of genomic locations. This is set to 1
        #if no value is passed. 
        if len(chromAccessibility) == 0:
            self.chromAccessibility = np.ones( self.N )
        else:
            self.chromAccessibility = chromAccessibility

        if len(secondTFintEnergies) > 0:
            locRange = np.arange(self.N)
            #Locations where the interaction energy is negative are
            #cooperatively bound by A and B
            coopLocations = locRange[ secondTFintEnergies < 0 ]

            #Locations where the interaction energy is positive are
            #competitively bound  by A and B.
            compLocations = locRange[ secondTFintEnergies > 0 ]

            #Binding energies of the second TF B. 
            self.locations.loc[:,'energy_B'] = secondTFspEnergies

            #Interaction energies between A and B at each genomic location.
            self.locations.loc[:,'int_energy'] = secondTFintEnergies

            #The binding type at each location is set to "cooperative" or "competitive"
            #at each location based on the interaction energy assigned to it. 
            self.locations.loc[coopLocations,'binding'] = 'cooperative'
            self.locations.loc[compLocations,'binding'] = 'competitive'

        if len( indirectLocations ) > 0:
            self.locations.loc[indirectLocations,'binding'] = 'indirect'
            self.indirectLocations = indirectLocations

        #pTFbound are the probabilities of finding the target TF bound at each
        #genomic location.  pBgBound is the probability of each of these
        #locations being bound in the input sample of the ChIP-seq experiment.
        pTFbound, pBgBound = self.computeBindingProbabilities( )
        self.locations.loc[:,'p_occ_bg'] = pBgBound

        chipFragments = binom.rvs( self.chipCells, self.locations['p_occ_chip'].values, size=self.N )
        self.locations.loc[:,'chip_fragments'] = chipFragments

        controlBound = binom.rvs( self.controlCells, self.locations['p_occ_bg'].values, size=self.N )
        self.locations.loc[:,'control_fragments'] = controlBound

    def computeBindingProbabilities( self ):
        """
        Compute binding probabilities at each genomic location based on its
        binding and interaction energy. 

        See Methods section in the manuscript for details on the calculation. 
        """

        unboundWt = np.exp( -self.unboundEnergy )

        #The probability of a location being bound in the input sample.
        bgWt = np.exp( -self.bgEnergy )
        pBgBound = bgWt/(unboundWt + bgWt)

        #The probability of a location being bound in the ChIP sample. 
        #This is the expression employed when there is only a single TF
        #capable of binding a location.
        spWt = np.exp( (-self.spEnergies + self.chemicalPotential) )
        pTFbound = spWt/(spWt + unboundWt)

        if len( self.secondTFintEnergies ) > 0 and len( self.indirectLocations ) == 0:
            #When there are two TFs present in the simulation, then the pTFbound
            #values computed in the earlier step are over-written. 
            #See the Methods section in the manuscript for the justification
            #behind computing occupancies of cooperatively bound TFs in this fashion.
            spWt = np.exp( (-self.spEnergies + self.chemicalPotential) )
            secondTFwt = np.exp( -self.secondTFspEnergies + self.secondTFchemicalPotential )
            coopWt = np.exp( -self.secondTFspEnergies - self.secondTFintEnergies - self.spEnergies + self.chemicalPotential + self.secondTFchemicalPotential )
            denom = (unboundWt + spWt + secondTFwt + coopWt)

            pTFbound = (spWt + coopWt)/denom
            
        if len( self.indirectLocations ) > 0:
            #In the case of indirect binding, the binding energy of the second
            #TF determines the occupancy of the location and not the binding
            #energy of the target TF.
            indirectWt = np.exp( -self.secondTFspEnergies + self.secondTFchemicalPotential )
            pTFbound[ self.indirectLocations ] = indirectWt[ self.indirectLocations ]/(unboundWt + indirectWt[self.indirectLocations])

        self.locations.loc[:,'p_occ_chip'] = pTFbound * self.chromAccessibility
        return [pTFbound,pBgBound]

