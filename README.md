# ANIMATE
A Python3 framework to simulate read counts in a ChIP-seq experiment. 

## Getting Started

Animate requires the following Python3 libraries installed --- ```argparse```,```numpy```,```scipy``` and ```pandas```. Animate has been tested to work with these versions of these libraries ---
    argparse : 1.1
    numpy : 1.15.2
    scipy : 1.1.0
    pandas : 0.23.4

These packages can be installed using the Python3 ```pip``` installer with the command ```pip3 install <package>```. 

## Running ANIMATE
	usage: animate.py [-h] [--mu-A MU_A] [--mu-B MU_B] [-c CONTROL_CELL_FRACTION]
                  [-b INPUT_BG] [-n NUM_CELLS] [-d DEPTH] [-p PCR_CYCLES] -i
                  INPUT_FILE [-o OUTPUT_PREFIX]

	The ANIMATE pipeline for simulating read counts in a ChIP-seq experiment

	Optional arguments:
  	-h, --help            show this help message and exit
  	--mu-A MU_A           Chemical potential (in units of k_B T) of TF A, where
        	                A is the target TF of the ChIP-seq. (default: 3.0)
  	--mu-B MU_B           Chemical potential (in units of k_B T) of TF B, where
                        B is a second TF that may be involved in cooperative
                        or indirect interactions with A, the target TF of the
                        ChIP-seq. (default: 3.0)
                        ChIP-seq. (default: 3.0)
	  -c CONTROL_CELL_FRACTION, --control-cell-fraction CONTROL_CELL_FRACTION
                        Control cell ratio. This is the fraction ofthe number
                        of cells used in the ChIP sample that is used for the
                        control sample. This value should be between 0 and 1.
                        Setting this parameter to 1 sets the number of cells
                        used in the ChIP and control samples to 1. (default:
                        0.1)
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
  	-o OUTPUT_PREFIX, --output-prefix OUTPUT_PREFIX
                        Name of the output file. The output is a tab separated
                        file that lists the following columns --- <chip_reads>
                        <unique_chip_reads> <control_reads>
                        <unique_control_reads>. See README for more
                        information on each column. In addition to the output 
			file above, two more output files are generated with the 
			suffixes '.run_info' and '.diag_output'. The run_info file 
			contains the parameters with which the output was generated 
			while the '.diag_output' contains results from the intermediate 
			steps of computation when the read counts are generated. (default: None)

	
## Input file
	
The header of the GENOME_FILE input must contain the following column items, of which ```p_ext```, ```p_amp``` and ```energy_A``` are the minimum columns that must be specified to run Animate (see the Examples section of the README). The GENOME_FILE columns represent the following quantities:

	p_ext --- The extraction efficiency at each genomic location. The value must
	lie between 0 and 1.
	p_amp --- PCR efficiency at each genomic location, which must lie between 0
	and 1. The mean number of amplified fragments at a location is (1 + p)^n, where
	p is the PCR efficiency at the location and n is the number of PCR cycles. Note
	that the PCR efficiencies are truncated to two decimal places in order to speed
	up the process of computing the number of amplified fragments obtained at each
	genomic location.  
	energy_A --- Binding energies of TF A at each location, which is the
	target TF of the ChIP-seq. A binding energy of 0 represents the strongest
	binding site with positive values representing weaker binding. 
	energy_B --- Binding energies of TF B at each location. This is only
	in case a second TF is being simulated.
	binding_type --- When no second TF is supplied, every location is considered
	to be directly bound by the TF A i.e., the value at each location defaults to
	"direct". If a second TF is supplied, the binding_type can be either be
	"direct" or "indirect". 
	int_energy--- When binding energies for A and B are supplied, the
	interaction energy at each location determines whether it is cooperatively,
	competitively or independently bound. At a given location, a negative value
	represents a cooperative interaction, a positive value represents a competitive
	interaction and a value of zero represents no interaction i.e., independent
	binding. 
	sequence--- The sequence of the region being bound. This sequence can be of
	any length. 
	chrom_accessibility--- The chromatin accessibility at each location. This is
	a value that lies between 0 and 1.

## Output generated	

The main output file, with the extension ```.animate.out``` is a tsv file that includes the following columns, in addition to the columns specified in the input file ---

	chip_reads --- The number of reads in the ChIP sample at each genomic
	location. This includes both unique reads and duplicate reads arising from
	sequencing PCR duplicate fragments during library preparation.
	unique_chip_reads --- The number of unique reads in the ChIP sample at each
	genomic location.
	control_reads --- The number of reads in the control sample at each genomic
	location. This includes both unique reads and duplicate reads arising from
	sequencing PCR duplicate fragments during library preparation.
	unique_control_reads --- The number of unique reads in the control sample at
	each genomic location.
	
Two more files are generated in addition to the ```.animate.out``` file, which are suffixed with ```.animate.out.run_info``` and ```.animate.out.diag_output```. The run_info file conains the parameters with which the output was generated  while the '.diag_output' contains results from the intermediate steps of computation when the read counts are generated. See
the Examples section for more details.

## Running the examples

### Minimum working example

The file ```basicExample.tsv``` in the ```examples``` folder is the minimum working input required to run animate.py. The contents of the file are the following ---

	p_ext	p_amp	energy_A
	0.539179	0.18	0.15
	0.505944	0.58	0.15
	0.498672	0.58	0.41
	0.479857	0.79	0.15
	0.494356	0.58	0.15
	0.554812	0.58	0.15
	0.545281	0.38	0.41
	0.492197	0.58	0.41
	0.503651	0.58	0.41
	0.578344	0.28	0.15

To execute ```animate.py``` on this file with the remaining setting set to their default values, run the following command ---

	python3 animate.py -i basicExample.tsv
	
No output file is specified in the above syntax. When this is the case, the input file name is used as a prefix to generate three output files. In this case, the files are named ```basicExample.tsv.animate.out``` (the main output file), ```basicExample.tsv.animate.out.diag_output``` (contains intermediate output from Animate along with the output from the main output file) and ```basicExample.tsv.animate.run_info``` (contains run information about the parameters of the simulation). 

In a run of the command ```python3 animate.py -i basicExample.tsv```, the contents of ```basicExample.tsv.animate.out``` are ---

	p_ext	p_amp	energy_A	chip_reads	unique_chip_reads	control_reads	unique_control_reads
	0.539179	0.18	0.15	0	0	0	0
	0.505944	0.58	0.15	84	70	82	66
	0.498672	0.58	0.41	84	66	98	81
	0.479857	0.79	0.15	481	185	487	187
	0.494356	0.58	0.15	73	57	76	65
	0.554812	0.58	0.15	71	58	77	69
	0.545281	0.38	0.41	11	11	6	6
	0.49219700000000005	0.58	0.41	81	66	75	64
	0.5036510000000001	0.58	0.41	109	80	93	71
	0.578344	0.28	0.15	6	5	6	6

The information on the parameters used to generate this output are in ```basicExample.tsv.animate.out.run_info``` ---

	Number of cells in ChIP sample : 100000
	Control cell ratio : 0.1
	Number of cells in control sample : 10000
	Chemical Potential of A : 3.0
	Number of PCR cycles : 15
	Sequencing depth : 100
	Total read count : 1000
	
The diagnostic output generated is in ```basicExample.tsv.animate.out.diag_output``` ---

	name	energy_A	binding	sequence	p_occ_chip	p_occ_bg	chip_fragments	control_fragments	unique_control_reads	control_reads	unique_chip_reads	chip_reads	amp_control_fragments	amp_chip_fragments	ext_control_fragments	ext_chip_fragments	read_count_ratio
	1	0.15	direct		0.9453186827840592	0.04742587317756678	94473	497	0	0	0	0	3396	2674	284	239	
	2	0.15	direct		0.9453186827840592	0.04742587317756678	94503	478	66	82	70	84	234729	204422	245	222	1.0606060606060606
	3	0.41	direct		0.9302152171234199	0.04742587317756678	93037	463	81	98	66	84	238212	227874	243	233	0.8148148148148148
	4	0.15	direct		0.9453186827840592	0.04742587317756678	94359	455	187	487	185	481	1298190	1352338	211	217	0.9893048128342246
	5	0.15	direct		0.9453186827840592	0.04742587317756678	94461	465	65	76	57	73	215390	214698	235	224	0.8769230769230769
	6	0.15	direct		0.9453186827840592	0.04742587317756678	94465	461	69	77	58	71	239417	251417	257	240	0.8405797101449275
	7	0.41	direct		0.9302152171234199	0.04742587317756678	93046	472	6	6	11	11	31739	35954	248	265	1.8333333333333333
	8	0.41	direct		0.9302152171234199	0.04742587317756678	93064	456	64	75	66	81	206097	221844	220	227	1.03125
	9	0.41	direct		0.9302152171234199	0.04742587317756678	93089	456	71	93	80	109	209127	273018	221	260	1.1267605633802817
	10	0.15	direct		0.9453186827840592	0.04742587317756678	94540	500	6	6	5	6	11240	11229	279	285	0.8333333333333334

### Running Animate with command-line parameters different from the default values

The sequencing depth, chemical potential, number of cells used in the control sample, and the number of PCR cycles can be changed from the command line. The remaining experimental parameters such as extraction efficiency, amplification efficiency need to be changed in the input file. 

The following command runs the file ```basicExample.tsv``` with different parameters and the output prefix set to a user-defined value ---

	python3 animate.py -i examples/basicExample.tsv --mu-A 1.5 --depth 300 --num-cells 10000 --control-cell-fraction 0.4 -o examples/mwe
	
The main output is now written to ```examples/mwe.animate.out``` ---

	p_ext	p_amp	energy_A	chip_reads	unique_chip_reads	control_reads	unique_control_reads
	0.539179	0.18	0.15	2	2	1	1
	0.505944	0.58	0.15	207	64	240	84
	0.498672	0.58	0.41	301	82	247	83
	0.479857	0.79	0.15	1368	69	1513	91
	0.494356	0.58	0.15	221	71	211	65
	0.554812	0.58	0.15	264	85	234	80
	0.545281	0.38	0.41	40	28	39	31
	0.49219700000000005	0.58	0.41	261	84	240	84
	0.5036510000000001	0.58	0.41	329	91	260	79
	0.578344	0.28	0.15	7	7	15	13
	
The file ```examples/mwe.animate.run_info``` shows the modified parameters ---

	Number of cells in ChIP sample : 10000
	Control cell ratio : 0.4
	Number of cells in control sample : 4000
	Chemical Potential of A : 1.5
	Number of PCR cycles : 15
	Sequencing depth : 300.0
	Total read count : 3000.0
	
The diagnostic output is written to ```examples/mwe.animate.diag_output``` ---

	name	energy_A	binding	sequence	p_occ_chip	p_occ_bg	chip_fragments	control_fragments	unique_control_reads	control_reads	unique_chip_reads	chip_reads	amp_control_fragments	amp_chip_fragments	ext_control_fragments	ext_chip_fragments	read_count_ratio
	1	0.15	direct		0.7941296281990528	0.04742587317756678	7964	169	1	1	2	2	1096	1431	90	102	2.0
	2	0.15	direct		0.7941296281990528	0.04742587317756678	7959	172	84	240	64	207	95250	65008	98	77	0.7619047619047619
	3	0.41	direct		0.7483817216070642	0.04742587317756678	7487	181	83	247	82	301	88268	84959	93	94	0.9879518072289156
	4	0.15	direct		0.7941296281990528	0.04742587317756678	7983	184	91	1513	69	1368	546708	436389	92	69	0.7582417582417582
	5	0.15	direct		0.7941296281990528	0.04742587317756678	7907	163	65	211	71	221	64847	73790	71	79	1.0923076923076922
	6	0.15	direct		0.7941296281990528	0.04742587317756678	7980	184	80	234	85	264	88031	90382	102	101	1.0625
	7	0.41	direct		0.7483817216070642	0.04742587317756678	7451	193	31	39	28	40	13153	10078	100	82	0.9032258064516129
	8	0.41	direct		0.7483817216070642	0.04742587317756678	7478	199	84	240	84	261	97850	88916	98	95	1.0
	9	0.41	direct		0.7483817216070642	0.04742587317756678	7551	187	79	260	91	329	88385	103121	93	106	1.1518987341772151
	10	0.15	direct		0.7941296281990528	0.04742587317756678	7963	191	13	15	7	7	4572	3334	102	84	0.5384615384615384

