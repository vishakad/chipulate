import numpy as np
from os import path
from itertools import product
from pandas import read_table

#This is the Tye7 energy matrix from the BEEML database, extracted from
#gr09.v2_Tye7-11_deBruijn.pwm.txt in the collection downloaded from
#http://stormo.wustl.edu/TF-BEMs/Mono-Nuc-Models.tar.gz 
#For more details on how #these matrices were derived, see the following
#publication --- Zhao, Yue, and Gary D. Stormo. "Quantitative analysis
#demonstrates most transcription factors require only simple models of
#specificity." Nature biotechnology 29.6 (2011): 480.
Tye7 = np.array( [[0.  , 1.39, 2.02, 0.  , 3.06, 1.59, 2.64, 7.96, 0.  , 2.46],
       [0.72, 0.85, 0.  , 2.46, 0.  , 1.89, 8.21, 7.23, 2.11, 0.7 ],
       [0.15, 1.27, 1.19, 0.84, 0.41, 0.  , 3.29, 0.  , 0.6 , 1.45],
       [1.64, 0.  , 1.84, 0.66, 1.88, 2.44, 0.  , 2.85, 1.45, 0.  ]])

#This is the Cbf1 energy matrix from the BEEML database, extracted from
#nbt06.v2_Cbf1_deBruijn_v2.pwm.txt in the collection downloaded from
#http://stormo.wustl.edu/TF-BEMs/Mono-Nuc-Models.tar.gz 
Cbf1 = np.array( [[0.42, 0.38, 1.87, 5.7, 0.0, 6.62, 2.31, 2.93, 4.59, 0.3],
[0.71, 1.79, 2.1, 0.0, 3.95, 0.0, 3.54, 3.06, 4.66, 0.0],
[0.0, 0.0, 1.81, 4.84, 5.75, 3.39, 0.0, 2.94, 0.0, 0.36],
[0.87, 2.0, 0.0, 6.6, 4.77, 3.13, 6.75, 0.0, 2.76, 0.14]] )

class MotifTable:
    def __init__( self, tf, pwm=[] ):
        """
        The MotifTable class encapsulates the energy matrix of a TF and provides
        routines for generating binding site sequences that fit a given binding
        energy distribution. 
        """
        if len(pwm) == 0 and tf.upper() not in ['TYE7','CBF1']:
            pwm = loadMotif( tf )
        elif tf.upper() == 'TYE7':
            pwm = Tye7 
        elif tf.upper() == 'CBF1':
            pwm = Cbf1

        self.pwm = pwm
        self.tf = tf

        if pwm.shape[0] == 4:
            pwm = pwm.transpose()

        self.maxEnergy = np.max( self.pwm, axis=1 ).sum()
        self.minEnergy = np.min( self.pwm, axis=1 ).sum()

    def getInformation( self, sequences ):
        """
        This function uses the PWM stored in the object and computes the motif
        scores of binding site sequences that are passed to it. 

        Input : 
        sequences --- An array of N sequences. Each sequence must be of length
        K, where K is the number of columns (or rows) of the binding energy
        matrix stored in the object.

        Returns an array of N values that are the energy/information of each of
        the sequences passed to the function.
        """
        locDict = {'A' : 0, 'C' : 1, 'G' : 2, 'T' : 3}
        if type(sequences) == 'str':
            sequences = [sequences]

        information = np.zeros( len(sequences) )

        if self.pwm.shape[0] == 4:
            pwm = self.pwm.transpose()
        else:
            pwm = self.pwm

        pwmLen = pwm.shape[0]
        pwmPositions = range( pwmLen )

        idx = 0
        for sequence in sequences:
            if 'N' in sequence:
                information[idx] = np.nan
                idx += 1
                continue

            convertedSeq = np.zeros( len(sequence), dtype=np.int )

            pos = 0
            for letter in sequence:
                convertedSeq[pos] = locDict[letter]
                pos += 1

            letterList = convertedSeq[0:pwmLen]
            information[idx] = pwm[pwmPositions,letterList].sum()
                
            idx += 1

        return information

    def sampleFromDistribution( self, energies ):
        """
        This function finds sequences whose binding energies closely match those
        of the energies passed to the function. 

        Input : 
        energies --- An array of N values that are positive mismatch binding
        energies. 
        
        Output : 
        There are two items returned as a list in the output.    
        The first item is an array of the exact energies of the sequences that
        were found.
        The second item is an array of the sequences found that matched the
        energies that were passed as an input.

        This function also writes to disk the binding energies of all possible binding sites
        as computed from the stored energy matrix.
        """
        revLocDict = {0 : 'A', 1 : 'C', 2 : 'G', 3 : 'T'}
        revLocArr = np.array(['A','C','G','T'])
        bsLength = self.pwm.shape[0]

        #outputEnergiesFile contains the energies of all possible binding sites,
        #as computed from the stored energy matrix.
        #outputSequencesFile contains the sequences of all possible binding sites.
        outputEnergiesFile = path.join( 'data', 'energies', '{}.all_energies.npy'.format( self.tf.upper() ) )
        outputSequencesFile = path.join( 'data', 'energies', '{}.all_sequences.npy'.format( self.tf.upper() ) )

        if self.pwm.shape[0] == 4:
            pwm = self.pwm.transpose()
        else:
            pwm = self.pwm
            
        pwmLen = pwm.shape[0]

        if path.isfile( outputEnergiesFile ) and path.isfile( outputSequencesFile ):
            allEnergies = np.load( outputEnergiesFile )
            allSequences = np.load( outputSequencesFile )
        else:
            alphabet = [0,1,2,3]
            pwmRange = range(pwmLen)
            #The product function from itertools helps generate all possible
            #binding site sequences
            sequenceIterator = product( alphabet, repeat=pwmLen)

            numSequences = 4 ** pwmLen
            allEnergies = np.zeros( numSequences, dtype=np.float )
            allSequences = np.zeros( (numSequences,pwmLen), dtype=np.uint8 )
            idx = 0
            for sequence in sequenceIterator:
                allSequences[idx] = sequence
                allEnergies[idx] = pwm[pwmRange,sequence].sum()
                idx += 1

            #The outputEnergiesFile and outputSequencesFile are created anew
            #in case no such file has been created for the TF associated
            #with the matrix. 
            np.save( outputEnergiesFile, allEnergies, allow_pickle=False )
            np.save( outputSequencesFile, allSequences, allow_pickle=False )

        mask = (energies >= self.minEnergy ) & (energies <= self.maxEnergy)
        if np.sum(~mask) > 0:
            print("{} sequences whose energies lie outside the range [{},{}] cannot be generated".format( np.sum(~mask), self.minEnergy, self.maxEnergy ))
        energies = energies[mask]

        #Once all possible binding sites and their corresponding energies
        #are computed, we pick those binding sites are closest in energy
        #to the energies passed into this function. 
        N = len(energies) 
        sampledEnergies = np.zeros( N )
        sampledSequences = np.zeros( N, dtype=np.object )
        sortedIdxes = np.argsort( allEnergies )
        allEnergies = allEnergies[sortedIdxes]
        allSequences = allSequences[sortedIdxes]

        itr = 0
        for energy in energies:
            idx = np.searchsorted( allEnergies, energy )
            sampledEnergies[itr] = allEnergies[idx]
            sampledSequences[itr] = "".join(revLocArr[allSequences[idx]])
            itr += 1

        return [sampledEnergies,sampledSequences]

def makeDinucBackground( numSequences, seqLen=10  ):
    """
    This function returns sequences whose dinucleotide frequencies match those
    of the S. cerevisiae genome (S288C/R64-2-1). 

    Input : 
    numSequences --- Number of sequences to be returned.
    seqLen --- The length of each sequence to be returned. This is a keyword
    argument whose default value is 10. 

    Returns :
    An array of sequences, each of length seqLen. 
    """
    dinucOrder = np.array([ 'AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT'])
    dinucFrequencies = [1.084e-01, 5.250e-02, 5.816e-02, 9.024e-02, 6.452e-02, 3.883e-02, 2.923e-02, 5.816e-02, 6.211e-02, 3.729e-02, 3.883e-02, 5.250e-02, 7.428e-02, 6.211e-02, 6.452e-02, 1.084e-01] 

    dinucTransMat = np.reshape( np.array( dinucFrequencies ), (4,4) )
    norm = np.repeat( dinucTransMat.sum(axis=1), 4 ).reshape(4,4)
    dinucTransMat = dinucTransMat/norm
    rowCumSums = np.cumsum( dinucTransMat, axis=1 )
    randVals = np.random.random( size=numSequences*seqLen )
    letters = np.array(['A','C','G','T'])

    sequences = np.zeros( numSequences*seqLen, dtype=np.object )

    firstLetterIdx = np.random.randint(4)
    firstLetter = letters[firstLetterIdx]
    sequences[0] = firstLetter
    prevLetter = firstLetter
    prevLetterIdx = firstLetterIdx
    for letterIdx in range(1,numSequences*seqLen):
        colToSelect = np.searchsorted( rowCumSums[prevLetterIdx,:], randVals[letterIdx]  )
        sequences[letterIdx] = letters[colToSelect]
        prevLetterIdx = colToSelect

    sequences = np.reshape( sequences, (numSequences,seqLen) )
    sequencesMerged = []
    for seq in sequences:
        sequencesMerged.append( "".join( seq ) )

    return np.array( sequencesMerged )

def klDist( reference, fromSim ):
    """
    This function computes the K-L distance between the two weight matrices
    that are passed to it.
    """
    kl = 0
    for (refRow,simRow) in zip(reference.transpose(),fromSim.transpose()):
        kl += -np.sum( refRow * np.log2( simRow/refRow ) )

    return kl

def loadMotif( tf ):
    #tf-filename.map contains the file names of energy matrices corresponding to
    #each TF passed to MotifTable. These are the matrices that have the highest
    #R^2 values with the microarray probe intensities from which they were inferred.
    df = read_table( path.join( 'data', 'tf-filename.map' ), sep="\t", header=None, names=['tf','filename'])
    row = df.query( 'tf == "{}"'.format( tf.upper() ) )
    if len(row) > 0:
        pwmFileName = row['filename'].values[0]
        df = read_table( path.join( 'data', 'pbm_pwms', pwmFileName ) + '.pwm.txt', sep=' ', header=None  )
        mat = df.loc[:,1:].values  #This gets rid of the first column, which is a positional index
        minEnergyMat = np.repeat( mat.min(axis=1), 4 ).reshape( (10,4) )
        mat -= minEnergyMat
        mat = mat.transpose()
        return mat
