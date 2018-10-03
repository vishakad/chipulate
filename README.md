# ANIMATE
A Python3 framework to simulate read counts in a ChIP-seq experiment. 

## Getting Started

Animate requires the following Python3 libraries installed --- ```argparse```,```numpy```,```scipy``` and ```pandas```. Animate has been tested to work with these versions of these libraries ---
    argparse : 1.1
    numpy : 1.15.2
    scipy : 1.1.0
    pandas : 0.23.4

## Running ANIMATE
	Usage: animate.py [-h] [--mu-A MU_A] [--mu-B MU_B] [-i INPUT_BG]
                  [-n NUM_CELLS] [-d DEPTH] [-p PCR_CYCLES] -g GENOME_FILE -o
                  OUTPUT_FILE

	The ANIMATE pipeline for simulating read counts in a ChIP-seq experiment

	Optional arguments:
  	-h, --help            show this help message and exit
  	--mu-A MU_A           Chemical potential (in units of k_B T) of TF A, where
        	                A is the target TF of the ChIP-seq. (default: 3.0)
  	--mu-B MU_B           Chemical potential (in units of k_B T) of TF B, where
                        B is a second TF that may be involved in cooperative
                        or indirect interactions with A, the target TF of the
                        ChIP-seq. (default: 3.0)
  	-i INPUT_BG, --input-bg INPUT_BG
                        Background binding energy (in units of k_BT) in the
                        input sample of the ChIP-seq experiment. Must be
                        greater than zero. A higher value indicates weaker
                        binding in the input sample. (default: 3.0)
  	-n NUM_CELLS, --num-cells NUM_CELLS
                        Number of cells used in the ChIP sample of the
                        experiment. A progressive increase in this value slows
                        down the simulation. (default: 100000)
  	-d DEPTH, --depth DEPTH
                        Sequencing depth. We define this as the number of
                        reads expected per binding location if an equal number
                        of reads came from each location. The total number of
                        sequence reads used is the product of the sequencing
                        depth and the number of binding locations. A
                        fractional value can be passed if the total number of
                        reads is less than the number of binding locations.
                        The coverage is set to be equal in both ChIP and input
                        samples. (default: 100)
  	-p PCR_CYCLES, --pcr-cycles PCR_CYCLES
                        Number of cycles employed in the PCR amplification
                        step. (default: 15)
  	-i GENOME_FILE, --input-file INPUT_FILE
                        File name of a tab-separated file that contains
                        location-wise information about the genome being
                        simulated and the experimental parameters of the ChIP-
                        seq. The first line is ignored and can be used as a header. 
			Each subsequent line contains an entry of the form <p_ext>
                        <p_amp> <binding_energy_A> <|sequence|> <|binding_energy_B|>
                        <|binding_type|> <|interaction energy|>
                        <|chrom_accessibility|>. The columns enclosed in
                        |..| are optional. See README for more information on
                        each column. (default: None)
  	-o OUTPUT_FILE, --output-file OUTPUT_FILE
                        Name of the output file. The output is a tab separated
                        file that lists the following columns --- <chip_reads>
                        <unique_chip_reads> <control_reads>
                        <unique_control_reads>. See README for more
                        information on each column. (default: None)

	The GENOME_FILE columns represent the following quantities:
	<p_ext> --- The extraction efficiency at each genomic location. The value must lie between 0 and 1.
	<p_amp> --- PCR efficiency at each genomic location, which must lie between 0 and 1. The mean number of amplified fragments at a location is (1 + p)^n, where p is the PCR efficiency at the location and n is the number of PCR cycles. Note that the PCR efficiencies are truncated to two decimal places in order to speed up the process of computing the number of amplified fragments obtained at each genomic location.  
	<binding_energy_A> --- Binding energies of TF A at each location, which is the target TF of the ChIP-seq. A binding energy of 0 represents the strongest binding site with positive values representing weaker binding. 
	<binding_energy_B> --- Binding energies of TF B at each location. This is only in case a second TF is being simulated.
	<binding_type> --- When no second TF is supplied, every location is considered to be directly bound by the TF A i.e., the value at each location defaults to "direct". If a second TF is supplied, the binding_type can be either be "direct" or "indirect". 
	<interaction_energy>--- When binding energies for A and B are supplied, the interaction energy at each location determines whether it is cooperatively, competitively or independently bound. At a given location, a negative value represents a cooperative interaction, a positive value represents a competitive interaction and a value of zero represents no interaction i.e., independent binding. 
	<sequence>--- The sequence of the region being bound. This sequence can be of any length. 
	<chrom_accessibility>--- The chromatin accessibility at each location. This is a value that lies between 0 and 1.

	The OUTPUT_FILE columns represent the following quantities:
	<chip_reads> --- The number of reads in the ChIP sample at each genomic location. This includes both unique reads and duplicate reads arising from sequencing PCR duplicate fragments during library preparation.
	<unique_chip_reads> --- The number of unique reads in the ChIP sample at each genomic location.
	<control_reads> --- The number of reads in the control sample at each genomic location. This includes both unique reads and duplicate reads arising from sequencing PCR duplicate fragments during library preparation.
	<unique_control_reads> --- The number of unique reads in the control sample at each genomic location.
