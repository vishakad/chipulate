#!/usr/bin/env python3
import os
import pandas as pd
import sys

import GenomeBindingTable as gbt
import FragExtract as Frag
import ChipSeq
import PCR

import numpy as np
import pandas as pd

import argparse
import pybedtools

def makeArray( val, N ):
    if type( val ) is list:
        val = np.array( val )
    elif type(val) not in [np.ndarray]:
        val = np.ones( N ) * val

    return val

def performChipSeq( sequences=[], spEnergies=[], numCells=100000, depth=100,
                   pExt=1.0, pAmp=0.58, pcrCycles=15, bgEnergy=1,
                   chemicalPotential=3, secondTFspEnergies=[],
                   secondTFchemicalPotential=0, chromAccessibility=[],
                   secondTFintEnergies=[], indirectLocations=[], controlCellRatio=1.0, generateIntervals=True ):
    """
    This function combines the GenomeBindingTable,FragExtract,PCR and ChIPseq
    classes together and returns a dataframe that contains read counts at
    each genomic location. 

    The inputs to this function are all keywords, which are listed below ---

    sequences --- Can be an empty list. Otherwise, a list of N sequences
    where N is the number of genomic locations. The sequences need not be
    equal in length. 

    spEnergies --- Binding energies of the target TF at each genomic
    location in the ChIP sample.

    bgEnergy --- Binding energy of each location in the input sample.
    chemicalPotential --- A value that is proportional the concentration
    the TF.  
    
    numCells --- Number of cells to be employed in the ChIP sample.

    depth --- This is the sequencing depth employed in the sample. The product
    of the number of locations and the sequencing depth determines the total
    number of reads.

    pExt --- The extraction efficiency of bound fragments in the ChIP and
    input samples. This can be either a scalar or an array of N value.

    pAmp --- The amplification efficiency of bound fragments in the ChIP and
    input samples. This can be either a scalar or an array of N value.

    pcrCycles --- The number of amplification cycles.

    controlCellRatio --- This is a fraction that determines the number of cells in the
    ChIP sample that will be employed in the control sample. The default
    value is 1.0 i.e. the same number of cells will be employed in both ChIP
    and input samples.

    secondTFspEnergies --- Binding energies of the second TF. 

    secondTFchemicalPotential --- Chemical potential of the second TF. 

    secondTFintEnergies --- An array that specifies the interaction
    energy between both TFs at each genomic location. Positive values
    indicate a competitive interaction, negative values indicate a
    cooperative interaction and zero indicates no interaction. 

    indirectLocations --- An array of location numbers that are
    to be simulated as being indirectly bound.  

    chromAccessibility -- An array of values that specify the chromatin
    accessibility at each genomic location. The values must lie between
    0 and 1. 

    Returns : 
    A dataframe with the following columns : 

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

    11) ext_chip_fragments --- The number of extracted fragments at each location
    in the ChIP sample

    12) ext_control_fragments --- The number of extracted fragments at each location
    in the control sample

    13) amp_chip_fragments --- The number of amplified fragments at each location
    in the ChIP sample

    14) amp_control_fragments --- The number of amplified fragments at each location
    in the control sample

    15) chip_reads --- The total number of reads at each genomic location in the
    ChIP sample.

    16) unique_chip_reads --- The total number of unique reads at each genomic
    location in the ChIP sample.

    15) control_reads --- The total number of reads at each genomic location in the
    control sample.

    16) unique_control_reads --- The total number of unique reads at each genomic
    location in the control sample.

    17) read_count_ratio --- The ratio between the number of unique ChIP reads
    and unique control reads at each genomic location. 
    """

    N = len( spEnergies )
    pExtChip = pExt
    pExtControl = pExt

    numChipReads = N * depth
    numControlReads = N * depth

    bgEnergy = makeArray( bgEnergy, N )
    pAmp = np.round( pAmp, 2 )
    pAmp = np.maximum( 0.01, pAmp )
    pAmp = np.minimum( pAmp, 0.99 )

    table = gbt.GenomeBindingTable( sequences, spEnergies,
                                   bgEnergy, chemicalPotential, numCells,
                                   secondTFspEnergies=secondTFspEnergies,
                                   secondTFchemicalPotential=secondTFchemicalPotential,
                                   secondTFintEnergies=secondTFintEnergies,
                                   indirectLocations=indirectLocations,
                                   controlCellRatio=controlCellRatio,
                                   chromAccessibility=chromAccessibility  )

    pExtControl = makeArray( pExtControl, N )
    pExtChip = makeArray( pExtChip, N )
    fragExtract = Frag.FragExtract( pExtControl, pExtChip, table )
    pAmp = makeArray( pAmp, N )
    pcrObj = ChipSeq.PCR( pcrCycles, pAmp )

    chipSeq = ChipSeq.ChipSeq( table, fragExtract, pcrObj, 
                              nControlReads=numControlReads,
                              nChipReads=numChipReads, generateIntervals=generateIntervals )
    
    genome = table.locations.merge( chipSeq.readsTable )
    genome = genome.merge( chipSeq.amplifiedTable )
    genome = genome.merge( fragExtract.extractedTable )
    genome.loc[:,'read_count_ratio'] = genome.eval('unique_chip_reads/unique_control_reads')

    return [genome,chipSeq.chipFragmentNumbers,chipSeq.controlFragmentNumbers]

def makeBed( bedDf, genome, chipFragmentNumbers, controlFragmentNumbers, chromSizesDf, fragmentLength=200, readLength=42, fragmentJitter=40, outputDir="", libraryType='single-end' ):
    fileNames = []
    outputPrefix = ""
    if len( outputDir ) == 0:
        outputPrefix = ""
    else:
        outputPrefix = outputDir + '.'

    for (fragmentStr,readsToDuplicate) in zip(['chip_reads','control_reads'],[chipFragmentNumbers,controlFragmentNumbers]):
        numFragments = genome[fragmentStr].sum()
        uniqueFragmentStr = 'unique_' + fragmentStr
        numUniqueFragments = genome[uniqueFragmentStr].sum()

        readsDf = pd.DataFrame( )
        readsDf.loc[:,'chr'] = np.zeros( numFragments, dtype=np.object )
        readsDf.loc[:,'start'] = np.zeros( numFragments, dtype=np.int64 )
        readsDf.loc[:,'end'] = np.zeros( numFragments, dtype=np.int64 )
        readsDf.loc[:,'name'] = np.zeros( numFragments, dtype=np.object )
        readsDf.loc[:,'score'] = np.ones( numFragments, dtype=np.int )
        readsDf.loc[:,'strand'] = np.zeros( numFragments, dtype=np.object )

        if libraryType == 'single-end':
            uniqueStrand = np.zeros( numUniqueFragments, dtype=np.object )
            uniqueStrandRand = np.random.random( size=numUniqueFragments )
            posStrand = uniqueStrandRand >= 0.5
            negStrand = ~posStrand
            uniqueStrand[posStrand] = '+'
            uniqueStrand[negStrand] = '-'

        regionNames = bedDf['name'].values
        if fragmentStr == 'chip_reads':
            uniqueFragmentMidPoints = np.repeat( bedDf['start'].values + bedDf['summit'].values, genome[uniqueFragmentStr])
            uniqueFragmentStartPos = np.random.normal( uniqueFragmentMidPoints - fragmentLength/2, fragmentJitter, size=numUniqueFragments ).astype(np.int64)
        elif fragmentStr == 'control_reads':
            uniqueFragmentMidPoints = np.repeat( bedDf['start'].values, genome[uniqueFragmentStr] )
            lengths = bedDf['end'].values - bedDf['start'].values
            lengths = np.repeat( lengths, genome[uniqueFragmentStr] )
            uniqueFragmentStartPos = np.int64( uniqueFragmentMidPoints + lengths * np.random.random( size=numUniqueFragments ) - fragmentLength/2 )
        
        uniqueReadNames = np.zeros( numUniqueFragments, dtype=np.object ) 

        idx = 0
        if libraryType == 'single-end':
            for (loc,count) in zip(range(bedDf.shape[0]),genome[uniqueFragmentStr]):
                uniqueReadNames[idx:(idx+count)] = [regionNames[loc] + '_read_{}'.format( num ) for num in np.arange(1,count+1)]
                idx += count
        elif libraryType == 'paired-end':
            uniqueReadNames2 = np.zeros( numUniqueFragments, dtype=np.object ) 
            for (loc,count) in zip(range(bedDf.shape[0]),genome[uniqueFragmentStr]):
                uniqueReadNames[idx:(idx+count)] = [regionNames[loc] + '_read_{}/1'.format( num ) for num in np.arange(1,count+1)]
                uniqueReadNames2[idx:(idx+count)] = [regionNames[loc] + '_read_{}/2'.format( num ) for num in np.arange(1,count+1)]
                idx += count
            
        fragmentStartPos = np.repeat( uniqueFragmentStartPos.astype(np.int64), readsToDuplicate )
        readNames = np.repeat( uniqueReadNames, readsToDuplicate )
        if libraryType == 'single-end':
            strand = np.repeat( uniqueStrand, readsToDuplicate )
            posStrand = strand == '+'
            negStrand = ~posStrand

        readsDf.loc[:,'chr'] = np.repeat( bedDf['chr'].values, genome[fragmentStr].values )
        readsDf = readsDf.merge( chromSizesDf, on='chr' )

        if libraryType == 'single-end':
            readsDf.loc[:,'strand'] = strand
            readsDf.loc[:,'name'] = readNames
            readsDf.loc[posStrand,'start'] = np.maximum( 1, fragmentStartPos[posStrand] )
            readsDf.loc[posStrand,'end'] = fragmentStartPos[posStrand] + readLength
            readsDf.loc[negStrand,'start'] = fragmentStartPos[negStrand] + fragmentLength - readLength
            readsDf.loc[negStrand,'end'] = fragmentStartPos[negStrand] + fragmentLength
        elif libraryType == 'paired-end':
            readsDf2 = readsDf.copy(deep=True)
            readNames2 = np.repeat( uniqueReadNames2, readsToDuplicate )

            readsDf.loc[:,'start'] = np.maximum( 1, fragmentStartPos )
            readsDf.loc[:,'end'] = fragmentStartPos + readLength 
            readsDf.loc[:,'strand'] = '+'
            readsDf.loc[:,'name'] = readNames

            readsDf2.loc[:,'start'] = fragmentStartPos + fragmentLength - readLength
            readsDf2.loc[:,'end'] = fragmentStartPos + fragmentLength
            readsDf2.loc[:,'strand'] = '-'
            readsDf2.loc[:,'name'] = readNames2

        outOfBounds = readsDf.query( 'end > max' ).index.tolist()
        readsLeftShift = 1 + np.random.randint( 10, size=len(outOfBounds) )
        readsDf.loc[outOfBounds,'end'] = readsDf.loc[outOfBounds,'max'] - readsLeftShift
        readsDf.loc[outOfBounds,'start'] = readsDf.loc[outOfBounds,'start'] - readsLeftShift
        if len( outOfBounds ) > 0:
            print( "For the following {} fragment(s) in the {} sample, their right ends were at positions beyond the length of the chromosome. They were shifted to the left by upto 10 base pairs. This can be avoided by either setting a shorter fragment length, or choosing genomic regions away from chromosome ends.".format( len(outOfBounds), fragmentStr[:-6] ) )
            print( readsDf[['chr','start','end','name','score','strand']].loc[outOfBounds,:] )

        if libraryType == 'single-end':
            fileName = outputPrefix + fragmentStr + '.bed'
            readsDf[['chr','start','end','name','score','strand']].to_csv( fileName,sep="\t",index=False,header=False)
            fileNames.append( fileName )
        else:
            for (ext,df) in zip(['_R1.bed','_R2.bed'],[readsDf,readsDf2]):
                fileName = outputPrefix + fragmentStr + ext
                print( fileName )
                df[['chr','start','end','name','score','strand']].to_csv( fileName,sep="\t",index=False,header=False)
                fileNames.append( fileName )

    return fileNames

def makeFastq( bedFileNames, genomeFileName, readLength, outputDir="", readErrors='none', libraryType='single-end' ):
    for bedFileName in bedFileNames:
        fastqFileName = bedFileName[:-4] + '.fastq'
        regionsFile = pybedtools.BedTool( bedFileName ).getfasta( fi=genomeFileName, bed=bedFileName, s=True, name=True )

        if outputDir == "":
            #regionsFile.save_seqs( fastaFileName )
            fastqFile = open( fastqFileName, 'w' )
        else:
            #regionsFile.save_seqs( '{}.{}'.format( outputDir, fastaFileName ) )
            fastqFile = open( '{}.{}'.format( outputDir, fastqFileName ), 'w' )


        fastaFile = open(regionsFile.seqfn)
        asciiBase = 33
        if readErrors == 'none':
            qualityScore = 42
            qualityString = chr( asciiBase + qualityScore ) * readLength
        
        for line in fastaFile:
            if line[0] == '>':
                line = line.replace('>','@')
                if libraryType == 'paired-end':
                    line = line[:-4] + '\n'
                fastqFile.write( line )
            else:
                fastqFile.write( line + '+\n{}\n'.format( qualityString ) )

        fastaFile.close()
        fastqFile.close()

parser = argparse.ArgumentParser(description='The ChIPulate pipeline for\
                                 simulating read counts in a ChIP-seq experiment', 
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter )

muAdefault = 3.0
parser.add_argument( '--mu-A', help='Chemical potential (in units of\
                             k_B T) of TF A, where A is the target TF of the\
                             ChIP-seq.'.format( muAdefault ), type=float, default=muAdefault)

muBdefault = 3.0
parser.add_argument( '--mu-B', help='Chemical potential (in units of\
                             k_B T) of TF B, where B is a second TF that may be\
                             involved in cooperative or indirect interactions\
                    with A, the target TF  of the ChIP-seq.'.format(
                            muBdefault), type=float, default=muBdefault,
                    required=False )

controlCellRatioDefault = 0.1
parser.add_argument( '-c', '--control-cell-fraction', help='Control cell ratio. This is the fraction of\
the number of cells used in the ChIP sample that is used for the control\
                    sample.  This value should be between 0 and 1. Setting this parameter\
                    to 1 sets the number of cells used in the ChIP and control\
                    samples to 1.', type=float, default=controlCellRatioDefault,
                    required=False )

bgEnergyDefault = 1.0
parser.add_argument( '-b', '--input-bg', help='Background binding energy (in units of\
                    k_BT) in the input sample of the ChIP-seq experiment. Must\
                    be less than the unbound energy. A higher value indicates weaker\
                    binding in the input sample.'.format( bgEnergyDefault ), type=float, default=bgEnergyDefault )

numCellsDefault = 100000
parser.add_argument( '-n', '--num-cells', help='Number of cells used in the ChIP\
                    sample of the experiment. A progressive increase in this value slows down\
                    the simulation.', type=int, default=numCellsDefault)

depthDefault = 100
parser.add_argument( '-d', '--depth', help='Sequencing depth. We define this as the\
                    number of reads expected per binding location if an equal\
                    number of reads came from each location. The total\
                    number of sequence reads used is the product of the\
                    sequencing depth and the number of binding locations. A\
                    fractional value can be passed if the total number of reads is\
                    less than the number of binding locations. The depth is\
                    set to be equal in both ChIP and input samples.', type=float, default=depthDefault )  

pcrCyclesDefault = 15
parser.add_argument( '-p', '--pcr-cycles', help='Number of cycles employed in the\
                    PCR amplification step.', type=int, default=pcrCyclesDefault )

parser.add_argument( '-i', '--input-file', help='File name of a tab-separated file that contains\
                    location-wise information about the genome being simulated\
                    and the experimental parameters of the ChIP-seq.\
                    Each line contains an entry of the form \
                    <p_ext> <p_amp> <binding_energy_A>    <|binding_energy_B|> <|binding_type|> <|interaction energy|> <|sequence|> <|chrom_accessibility|>, where \
                    the columns enclosed in |..| are optional. See README for more information on each column.', type=str, required=True )

chromSizesFileName = ""
parser.add_argument( '--chrom-size-file', help='File containing sizes of each chromosome. This argument is ignored when the chr, start and end columns are not specified in the input file. If these columns are specified, a tab-separated file where the first two columns contain the chromosome name and size, respectively, must be supplied. ', type=str, required=False, default=chromSizesFileName )

genomeFileName = ""
parser.add_argument( '-g', '--genome-file', help='File containing a genome sequence in FASTA format. If FASTA output is requested (by specifying the chr, start and end columns in the input file), a single FASTA file containing the genome sequence must be specified. If chr, start and end columns are not specified in the input, then this argument is ignored.', type=str, required=False, default=genomeFileName )

parser.add_argument( '-o', '--output-prefix', help='Prefix of the output file. The\ output is a tab separated file that lists the following\ columns --- <chip_reads> <unique_chip_reads> <control_reads>\ <unique_control_reads>. See README for more information on\ each column.', type=str, required=False, default=None )

parser.add_argument( '--output-dir', help='Directory to which all output should be written. Ensure that you have write privileges to this directory.', type=str, required=False, default=None )

defaultReadLength = 150
parser.add_argument( '--read-length', help='Read length (in base pairs) to simulate. This must be smaller than the fragment length(s) specified in the --fragment-length argument, and is a required argument if FASTQ output is requested. Only single-end reads are simulated.', type=int, required=False, default=defaultReadLength)

defaultFragmentLength = 200
parser.add_argument( '--fragment-length', help='Fragment length (in base pairs) to simulate. This must be larger than the read length specified for --read-length.', type=int, required=False, default=defaultFragmentLength)

defaultFragmentJitter = 20
parser.add_argument( '--fragment-jitter', help='Variation in the starting position of fragments (in base pairs) at a genomic region. A larger value leads to a greater variation in start positions of fragments in the ChIP sample only. This parameter does not affect the starting position of fragments in the control sample.', type=int, required=False, default=defaultFragmentJitter)

defaultLibraryType = 'single-end'
parser.add_argument( '--library-type', help='Type of sequencing library. This can be set to either "single-end" or "paired-end". If single-end is specified, a single BED and FASTQ file is generated for ChIP and control samples. If "paired-end" is specified, paired-end reads are simulated. In this case, two sets of BED and FASTQ files are produced for each of the ChIP and control samples.The BED and FASTQ files containing reads from the forward and reverse strands are suffixed with _R1 and _R2, respectively.', required=False, default=defaultLibraryType )

args = parser.parse_args()

def validateInput( df, args ):
    numLocations = df.shape[0]
    terminateFlag = False
    allowedColumnNames = ['chr','start','end','name','summit','p_ext','p_amp','energy_A','energy_B','sequence','binding_type','int_energy']
    hasOnlyNaNs = []

    #Check column names
    for column in df.columns:
        if column not in allowedColumnNames:
            print("ERROR : Column {} is not an allowed column name. Perhaps the entries are not tab separated?".format( column ) ,file=sys.stderr)
            terminateFlag = True

    if terminateFlag:
        print("ERROR : Allowed column names are {}, where ['p_ext','p_amp','energy_A'] are required inputs.".format( allowedColumnNames ) ,file=sys.stderr)

    #Check if correct datatypes have been specified in each column. 
    for column in df.columns:
        isNumericType = pd.api.types.is_numeric_dtype( df[column] )
        if column in ['p_ext','p_amp','energy_A','energy_B','int_energy'] and not isNumericType:
            print("ERROR : Only numeric data expected in column {}.".format( column ),file=sys.stderr)
            terminateFlag = True
        elif column in ['binding_type','sequence'] and isNumericType:
            print("ERROR : Only non-numeric data expected in column {}.".format( column ) ,file=sys.stderr)
            terminateFlag = True

    #Check if there are missing entries in any of the columns.
    columnsWithNaNs = []
    for column in df.columns:
        numNaNs = df[column].isnull().sum()
        if not (numNaNs == numLocations or numNaNs == 0):
            print( "ERROR : Column {} has missing entries.".format( column ) ,file=sys.stderr)
            terminateFlag = True
        elif numNaNs == numLocations:
            if column in ['p_ext','p_amp','energy_A']:
                print("ERROR : Values for column {} needs to be specified as it is a required input.".format( column ) ,file=sys.stderr)
                terminateFlag = True
            columnsWithNaNs.append( column )

    if ('energy_B' in df.columns and 'int_energy' not in df.columns): 
        print('ERROR : Interaction energies need to be specified if a second TF is to be simulated',file=sys.stderr)
        terminateFlag = True

    df = df.drop( columnsWithNaNs, axis=1 )
    
    #Check if entries in the p_ext, p_amp and chrom_accessibility columns lie
    #between 0 and 1
    acceptableRanges = {'p_ext' : [0,1], 'p_amp' : [0,1], 'energy_A' : [0,np.inf], 'energy_B' : [0,np.inf], 'chrom_accessibility' : [0,1], 'binding_type' : ['direct','indirect'], 'sequence' : set( ['A','a','C','c','G','g','T','t','N'])}
    for column in df.columns:
        if column in acceptableRanges.keys():
            if column not in ['binding_type','sequence']:
                maxVal = df[column].max()
                minVal = df[column].min()
                if minVal < acceptableRanges[column][0] or maxVal > acceptableRanges[column][1]:
                    print("ERROR : Some value(s) in column {} lie outside the range {}".format( column, acceptableRanges[column] ) ,file=sys.stderr)
                    terminateFlag = True 
            elif column == 'binding_type':
                bindingTypeSum = 0
                for bindingType in acceptableRanges[column]:
                    bindingTypeSum += df.query( 'binding_type == "{}"'.format( bindingType ) ).shape[0]
                
                if bindingTypeSum != numLocations:
                    print("ERROR : Some value(s) in column {} are not in the acceptable set of value {}".format( column, acceptableRanges[column] ) ,file=sys.stderr)
                    terminateFlag = False
            elif column == 'sequence':
                allSequences = "".join( df['sequence'].values )
                lettersPresent = set( allSequences )
                numLettersPresent = len( lettersPresent )

                intersection = lettersPresent.intersection( acceptableRanges[column] )
                if intersection != lettersPresent: 
                    print("WARNING : Invalid letters present in input sequences.  Accepted letters are {}. Continuing anyway".format( acceptableRanges[column] ) ,file=sys.stderr)
                    terminateFlag = True

    if 'int_energy' in df.columns and 'binding_type' in df.columns:
        numIndirectInteractionClashes = df.query('binding_type == "indirect" & int_energy > 0').shape[0]
        if numIndirectInteractionClashes > 0:
            print("Note : There are {} entries where the target TF is specified to indirectly bind DNA but an interaction energy has also been specified. These entries will be treated as indirectly bound by A via B and the interaction energy will be ignored.".format( numIndirectInteractionClashes ))

    return [terminateFlag,df]

def validateBedFasta( args ):
    chromSizesFileName = args.chrom_size_file
    genomeFileName = args.genome_file
    readLength = args.read_length
    fragmentLength = args.fragment_length
    fragmentJitter = args.fragment_jitter
    libraryType = args.library_type

    terminateFlag = False

    if len(chromSizesFileName) == 0:
        print("A chrom.sizes file should be passed in order to generate BED/FASTA output.", file=sys.stderr,file=sys.stderr)
        terminateFlag = True
    elif not os.path.isfile(chromSizesFileName):
        print("The chrom.sizes file at {} could not be found.".format(chromSizesFileName), file=sys.stderr,file=sys.stderr)
        terminateFlag = True
    else:
        chromSizesDf = pd.read_table( chromSizesFileName, sep="\t", header=None )
        if not np.issubdtype( chromSizesDf[1].dtype, np.number ):
            print( "Second column of chrom.sizes file should contain numeric data.", file=sys.stderr,file=sys.stderr)
            terminateFlag = True
        elif chromSizesDf.shape[0] == 0:
            print( "Chrom.sizes file appears empty." ,file=sys.stderr)
            terminateFlag = True

    if len(genomeFileName) == 0:
        print("A genome file should be passed in order to generate FASTA output.", file=sys.stderr,file=sys.stderr)
        terminateFlag = True
    elif not os.path.isfile(genomeFileName):
        print("The genome file at {} could not be found.".format(genomeFileName),file=sys.stderr)
        terminateFlag = True

    if libraryType not in ['single-end','paired-end']:
        print("Library type was specified as {}. This should either be \"single-end\" or \"paired-end\".".format( libraryType ) ,file=sys.stderr)
        terminateFlag = True

    for (val,quantity) in zip( [fragmentJitter,fragmentLength,readLength], ['Fragment jitter','Fragment length','Read length'] ):
        if val < 0:
            print("{} was specified as {} bp. This should be a non-negative value.".format( quantity, val ),file=sys.stderr)
            terminateFlag = True

    if readLength > 0 and fragmentLength > 0 and fragmentLength < readLength:
        print("Fragment length specified ({} bp) is lower than the read length specified ({} bp). Read length must be less than fragment length.".format( fragmentLength, readLength), file=sys.stderr ,file=sys.stderr) 

    return terminateFlag 

def main():
    inputFileName = args.input_file
    if args.output_prefix is None:
        outputFileName = inputFileName + '.chipulate.out'
        outputPrefix = inputFileName
    else:
        outputFileName = args.output_prefix + '.chipulate.out'
        outputPrefix = args.output_prefix

    if args.output_dir is None:
        outputDir = outputPrefix
    else:
        outputDir = os.path.join( args.output_dir, outputPrefix )

    diagOutputFileName = outputDir + '.chipulate.diag_output'
    runInfoOutputFileName = outputDir + '.chipulate.run_info'

    #inputDf = pd.read_csv( inputFileName, sep="\t", skiprows=1, names=['p_ext','p_amp','energy_A','sequence','binding_type','energy_B','int_energy','chrom_accessibility'])
    inputDf = pd.read_csv( inputFileName, sep="\t" )
    numLocations = inputDf.shape[0]

    terminateFlagInput, inputDf = validateInput( inputDf, args ) 
    if 'chr' in inputDf.columns and 'start' in inputDf.columns and 'end' in inputDf.columns:
        generateIntervals = True
        terminateFlagBedFasta = validateBedFasta( args )
    
    depth = args.depth
    numCells = args.num_cells
    chemicalPotentialA = args.mu_A
    chemicalPotentialB = args.mu_B
    pcrCycles = args.pcr_cycles
    controlCellRatio = args.control_cell_fraction
    inputBgEnergy = args.input_bg
    chromSizesFileName = args.chrom_size_file
    genomeFileName = args.genome_file
    readLength = args.read_length
    fragmentLength = args.fragment_length
    fragmentJitter = args.fragment_jitter
    libraryType = args.library_type

    if terminateFlagBedFasta or terminateFlagInput :
        print("Error encountered in input. Aborting.", file=sys.stderr)
        return 0

    spEnergies = inputDf['energy_A']
    pAmp = inputDf['p_amp']
    pExt = inputDf['p_ext']
    if 'energy_B' in inputDf.columns:
        secondTFspEnergies = inputDf['energy_B']
        intEnergy = inputDf['int_energy']
    else:
        secondTFspEnergies = []
        secondTFintEnergies = []

    if 'binding_type' in inputDf.columns:
        allLocations = np.arange( numLocations )
        indirectLocations = allLocations[ inputDf['binding_type'] == 'indirect' ]
    else:
        indirectLocations = []

    if 'sequence' in inputDf.columns:
        sequences = inputDf['sequence']
    else:
        sequences = []

    if 'chrom_accessibility' in inputDf.columns:
        chromAccessibility = inputDf['chrom_accessibility']
    else:
        chromAccessibility = []

    bedCols = []
    if generateIntervals:
        bedCols = ['chr','start','end','name','summit','strand']

        #Assign random summits to each region.
        if 'summit' not in inputDf.columns:
            starts = inputDf['start'].values
            ends = inputDf['end'].values
            inputDf.loc[:,'summit'] = inputDf.eval( '(end-start)/2' )

        inputDf.loc[:,'strand'] = '.'
        chromSizesDf = pd.read_table( chromSizesFileName, sep="\t", header=None )
        chromSizesDf = chromSizesDf[[0,1]].rename( {0 : 'chr', 1 : 'max'}, axis=1 )

    outputDf, chipFragmentNumbers, controlFragmentNumbers = performChipSeq( sequences=sequences, spEnergies=spEnergies,
                            numCells=numCells, depth=depth, pAmp=pAmp,
                            pExt=pExt, pcrCycles=pcrCycles,
                            bgEnergy=inputBgEnergy, controlCellRatio=controlCellRatio,
                            chemicalPotential=chemicalPotentialA,
                            secondTFspEnergies=secondTFspEnergies,
                            secondTFchemicalPotential=chemicalPotentialB, 
                            chromAccessibility=chromAccessibility,
                            indirectLocations=indirectLocations, generateIntervals=generateIntervals )

    if 'name' not in inputDf.columns:
        inputDf.loc[:,'name'] = outputDf['name'].values

    if generateIntervals:
        bedFileNames = makeBed( inputDf[bedCols], outputDf, chipFragmentNumbers, controlFragmentNumbers, chromSizesDf, outputDir=outputDir, readLength=readLength, fragmentLength=fragmentLength, fragmentJitter=fragmentJitter, libraryType=libraryType )

        makeFastq( bedFileNames, genomeFileName, readLength, libraryType=libraryType )


    colsToWrite = inputDf.columns.tolist()
    for el in ['chip_reads','unique_chip_reads','control_reads','unique_control_reads']:
        colsToWrite.append( el )

    for col in outputDf.columns:
        inputDf.loc[:,col] = outputDf[col]

    #Write basic output to output file
    inputDf[colsToWrite].to_csv( outputFileName, sep="\t", index=False )

    #Write diagnostic output
    outputDf.to_csv( diagOutputFileName, sep="\t", index=False )

    #Write run information
    runInfoFile = open( runInfoOutputFileName, 'w' )

    runInfoFile.write( "Number of cells in ChIP sample : {}\n".format( numCells ) )
    runInfoFile.write( "Control cell ratio : {}\n".format( controlCellRatio ) )
    runInfoFile.write( "Number of cells in control sample : {}\n".format( np.int(controlCellRatio*numCells ) ) )
    runInfoFile.write( "Chemical Potential of A : {}\n".format( chemicalPotentialA) )
    if 'energy_B' in inputDf.columns:
        runInfoFile.write( "Chemical Potential of B : {}\n".format( chemicalPotentialB) )
    runInfoFile.write( "Number of PCR cycles : {}\n".format( pcrCycles ) )
    runInfoFile.write( "Sequencing depth : {}\n".format( depth ) )
    runInfoFile.write( "Total read count : {}\n".format( depth*numLocations ) )
    if generateIntervals:
        runInfoFile.write( "Library type : {}\nRead length : {}\nFragment length : {}\nFragment jitter : {}\n".format( args.library_type, args.read_length, args.fragment_length, args.fragment_jitter ) )

if __name__ == "__main__":
    main()
