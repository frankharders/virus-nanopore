#!/bin/bash

##  Shufflon analysis for nanopore data "fasta only!!!"
##  QC of input data 

##  activate the environment for this downstream analysis
eval "$(conda shell.bash hook)";
conda activate NGS-bbmap;

## create directories
mkdir -p ./00_fastqc;
mkdir -p ./LOGS;

## variables used
QC="$PWD"/00_fastqc;
RAW="$PWD"/RAWREADS;
LOG="$PWD"/LOGS;
NODES=48;

echo -e "\n\n\n de duplicate fastq file starting\n";

##  file with sample names
SAMPLEFILE="$PWD"/samples.txt;

#count0=1;
#countD=$(cat "$SAMPLEFILE" | wc -l);

#	while [ "$count0" -le "$countD" ];do

#		SAMPLE=$(cat "$SAMPLEFILE" | awk 'NR=='"$count0");

#echo -e "$SAMPLE";

# files
#		FILEin="$RAW"/"$SAMPLE".fastq.gz;
#		FILEout="$RAW"/"$SAMPLE".dedupe.fastq.gz;
#		LOG="$PWD"/"LOGS"/"$SAMPLE".general.log;

## deduplicate the original fastq files for duplicate reads
#			dedupe.sh in="$FILEin" out="$FILEout" qin=33 ignorejunk ow > "$LOG" 2>&1;

#		count0=$((count0+1));
#	done

#####

##  activate the environment for this downstream analysis
eval "$(conda shell.bash hook)";
conda activate NGS-porechop;

echo -e "\n\n\nadapter trimming is started\n";

##  loop started with sample names as input

count1=1;
countS=$(cat "$SAMPLEFILE" | wc -l);

	while [ "$count1" -le "$countS" ];do

		SAMPLE=$(cat "$SAMPLEFILE" | awk 'NR=='"$count1");

## files
		FILEin="$RAW"/"$SAMPLE".fastq.gz;
		FILEout="$RAW"/"$SAMPLE".chopped.fastq.gz;
		LOG="$PWD"/"LOGS"/"$SAMPLE".general.log;
##  adapter trimming
#			porechop_abi --ab_initio --format auto -t "$NODES" -i "$FILEin" -o "$FILEout" >> "$LOG" 2>&1 ;
			porechop -i "$FILEin" --format auto -o "$FILEout" -t "$NODES"  >> "$LOG" 2>&1;


	count1=$((count1+1));
	done

#####

##  activate the environment for this downstream analysis
eval "$(conda shell.bash hook)";
conda activate NGS-QC;

echo -e "\n\n\n fastqc of input files is started\n";

##  loop started with sample names as input


count2=1;
countS=$(cat "$SAMPLEFILE" | wc -l);

	while [ "$count2" -le "$countS" ];do

			SAMPLE=$(cat "$SAMPLEFILE" | awk 'NR=='"$count2");

## files
			FILEin="$RAW"/"$SAMPLE".chopped.fastq.gz;
			LOG="$PWD"/"LOGS"/"$SAMPLE".general.log;

##  QC for reporting
				fastqc -t "$NODES" --nano -o $QC $FILEin >> "$LOG" 2>&1;



		count2=$((count2+1));
	done


exit 1


##### porechop_abi #####
#
#usage: porechop_abi [-abi] [-go] [-abc AB_INITIO_CONFIG] [-tmp TEMP_DIR] [-cap CUSTOM_ADAPTERS] [-ddb] [-ws WINDOW_SIZE] [-ndc] [-nr NUMBER_OF_RUN] [-cr CONSENSUS_RUN]
#                    [-ec EXPORT_CONSENSUS] [-aax ALL_ABOVE_X] [-box BEST_OF_X] [--export_graph EXPORT_GRAPH] -i INPUT [-o OUTPUT] [--format {auto,fasta,fastq,fasta.gz,fastq.gz}]
#                    [-v VERBOSITY] [-t THREADS] [-b BARCODE_DIR] [--barcode_threshold BARCODE_THRESHOLD] [--barcode_diff BARCODE_DIFF] [--require_two_barcodes] [--untrimmed]
#                    [--discard_unassigned] [--adapter_threshold ADAPTER_THRESHOLD] [--check_reads CHECK_READS] [--scoring_scheme SCORING_SCHEME] [--end_size END_SIZE]
#                    [--min_trim_size MIN_TRIM_SIZE] [--extra_end_trim EXTRA_END_TRIM] [--end_threshold END_THRESHOLD] [--no_split] [--discard_middle]
#                    [--middle_threshold MIDDLE_THRESHOLD] [--extra_middle_trim_good_side EXTRA_MIDDLE_TRIM_GOOD_SIDE] [--extra_middle_trim_bad_side EXTRA_MIDDLE_TRIM_BAD_SIDE]
#                    [--min_split_read_size MIN_SPLIT_READ_SIZE] [-h] [--version]
#
#Porechop_ABI v_0.5.0: Ab Initio version of Porechop. A tool for finding adapters in Oxford Nanopore reads, trimming them from the ends and splitting reads with internal adapters
#
#Ab-Initio options:
#  -abi, --ab_initio                     Try to infer the adapters from the read set instead of just using the static database. (default: False)
#  -go, --guess_adapter_only             Just display the inferred adapters and quit. (default: False)
#  -abc AB_INITIO_CONFIG, --ab_initio_config AB_INITIO_CONFIG
#                                        Path to a custom config file for the ab_initio phase (default file in Porechop folder)
#  -tmp TEMP_DIR, --temp_dir TEMP_DIR    Path to a writable temporary directory. Directory will be created if it does not exists. Default is ./tmp
#  -cap CUSTOM_ADAPTERS, --custom_adapters CUSTOM_ADAPTERS
#                                        Path to a custom adapter text file, if you want to manually submit some.
#  -ddb, --discard_database              Ignore adapters from the Porechop database. This option require either ab-initio (-abi) or a custom adapter (-cap) to be set. (default: False)
#  -ws WINDOW_SIZE, --window_size WINDOW_SIZE
#                                        Size of the smoothing window used in the drop cut algorithm. (set to 1 to disable). (default: 3)
#  -ndc, --no_drop_cut                   Disable the drop cut step entirely (default: False)
#
#Consensus mode options:
#  -nr NUMBER_OF_RUN, --number_of_run NUMBER_OF_RUN
#                                        Number of time the core module must be run to generate the first consensus. Each count file is exported separately. Set to 1 for single run
#                                        mode. (default: 10)
#  -cr CONSENSUS_RUN, --consensus_run CONSENSUS_RUN
#                                        With -nr option higher than 1, set the numberof additional runs performed if no stable consensus is immediatly found. (default: 20)
#  -ec EXPORT_CONSENSUS, --export_consensus EXPORT_CONSENSUS
#                                        Path to export the intermediate adapters found in consensus mode.
#  -aax ALL_ABOVE_X, --all_above_x ALL_ABOVE_X
#                                        Only select consensus sequences if they are made using at least x percent of the total adapters. Default is 10%.
#  -box BEST_OF_X, --best_of_x BEST_OF_X
#                                        Only select the best x consensus sequences from all consensus found. (default: 0)
#
#Graphs options:
#  --export_graph EXPORT_GRAPH           Path to export the assembly graphs (.graphml format), if you want to keep them
#
#Main options:
#  -i INPUT, --input INPUT               FASTA/FASTQ of input reads or a directory which will be recursively searched for FASTQ files (required)
#  -o OUTPUT, --output OUTPUT            Filename for FASTA or FASTQ of trimmed reads (if not set, trimmed reads will be printed to stdout)
#  --format {auto,fasta,fastq,fasta.gz,fastq.gz}
#                                        Output format for the reads - if auto, the format will be chosen based on the output filename or the input read format (default: auto)
#  -v VERBOSITY, --verbosity VERBOSITY   Level of progress information: 0 = none, 1 = some, 2 = lots, 3 = full - output will go to stdout if reads are saved to a file and stderr if
#                                        reads are printed to stdout (default: 1)
#  -t THREADS, --threads THREADS         Number of threads to use for adapter alignment (default: 16)
#
#Barcode binning settings:
#  Control the binning of reads based on barcodes (i.e. barcode demultiplexing)
#
#  -b BARCODE_DIR, --barcode_dir BARCODE_DIR
#                                        Reads will be binned based on their barcode and saved to separate files in this directory (incompatible with --output)
#  --barcode_threshold BARCODE_THRESHOLD
#                                        A read must have at least this percent identity to a barcode to be binned (default: 75.0)
#  --barcode_diff BARCODE_DIFF           If the difference between a read's best barcode identity and its second-best barcode identity is less than this value, it will not be put in a
#                                        barcode bin (to exclude cases which are too close to call) (default: 5.0)
#  --require_two_barcodes                Reads will only be put in barcode bins if they have a strong match for the barcode on both their start and end (default: a read can be binned
#                                        with a match at its start or end)
#  --untrimmed                           Bin reads but do not trim them (default: trim the reads)
#  --discard_unassigned                  Discard unassigned reads (instead of creating a "none" bin) (default: False)
#
#Adapter search settings:
#  Control how the program determines which adapter sets are present
#
#  --adapter_threshold ADAPTER_THRESHOLD
#                                        An adapter set has to have at least this percent identity to be labelled as present and trimmed off (0 to 100) (default: 90.0)
#  --check_reads CHECK_READS             This many reads will be aligned to all possible adapters to determine which adapter sets are present (default: 10000)
#  --scoring_scheme SCORING_SCHEME       Comma-delimited string of alignment scores: match, mismatch, gap open, gap extend (default: 3,-6,-5,-2)
#
#End adapter settings:
#  Control the trimming of adapters from read ends
#
#  --end_size END_SIZE                   The number of base pairs at each end of the read which will be searched for adapter sequences (default: 150)
#  --min_trim_size MIN_TRIM_SIZE         Adapter alignments smaller than this will be ignored (default: 4)
#  --extra_end_trim EXTRA_END_TRIM       This many additional bases will be removed next to adapters found at the ends of reads (default: 2)
#  --end_threshold END_THRESHOLD         Adapters at the ends of reads must have at least this percent identity to be removed (0 to 100) (default: 75.0)
#
#Middle adapter settings:
#  Control the splitting of read from middle adapters
#
#  --no_split                            Skip splitting reads based on middle adapters (default: split reads when an adapter is found in the middle)
#  --discard_middle                      Reads with middle adapters will be discarded (default: reads with middle adapters are split) (required for reads to be used with Nanopolish,
#                                        this option is on by default when outputting reads into barcode bins)
#  --middle_threshold MIDDLE_THRESHOLD   Adapters in the middle of reads must have at least this percent identity to be found (0 to 100) (default: 90.0)
#  --extra_middle_trim_good_side EXTRA_MIDDLE_TRIM_GOOD_SIDE
#                                        This many additional bases will be removed next to middle adapters on their "good" side (default: 10)
#  --extra_middle_trim_bad_side EXTRA_MIDDLE_TRIM_BAD_SIDE
#                                        This many additional bases will be removed next to middle adapters on their "bad" side (default: 100)
#  --min_split_read_size MIN_SPLIT_READ_SIZE
#                                        Post-split read pieces smaller than this many base pairs will not be outputted (default: 1000)
#
#Help:
#  -h, --help                            Show this help message and exit
#  --version                             Show program's version number and exit
#
#####

##### dedupe #####
#
#Written by Brian Bushnell and Jonathan Rood
#Last modified February 19, 2020
#
#Description:  Accepts one or more files containing sets of sequences (reads or scaffolds).
#Removes duplicate sequences, which may be specified to be exact matches, subsequences, or sequences within some percent identity.
#Can also find overlapping sequences and group them into clusters.
#Please read bbmap/docs/guides/DedupeGuide.txt for more information.
#
#Usage:     dedupe.sh in=<file or stdin> out=<file or stdout>
#
#An example of running Dedupe for clustering short reads:
#dedupe.sh in=x.fq am=f ac=f fo c pc rnc=f mcs=4 mo=100 s=1 pto cc qin=33 csf=stats.txt pattern=cluster_%.fq dot=graph.dot
#
#Input may be fasta or fastq, compressed or uncompressed.
#Output may be stdout or a file.  With no output parameter, data will be written to stdout.
#If 'out=null', there will be no output, but statistics will still be printed.
#You can also use 'dedupe <infile> <outfile>' without the 'in=' and 'out='.
#
#I/O parameters:
#in=<file,file>        A single file or a comma-delimited list of files.
#out=<file>            Destination for all output contigs.
#pattern=<file>        Clusters will be written to individual files, where the '%' symbol in the pattern is replaced by cluster number.
#outd=<file>           Optional; removed duplicates will go here.
#csf=<file>            (clusterstatsfile) Write a list of cluster names and sizes.
#dot=<file>            (graph) Write a graph in dot format.  Requires 'fo' and 'pc' flags.
#threads=auto          (t) Set number of threads to use; default is number of logical processors.
#overwrite=t           (ow) Set to false to force the program to abort rather than overwrite an existing file.
#showspeed=t           (ss) Set to 'f' to suppress display of processing speed.
#minscaf=0             (ms) Ignore contigs/scaffolds shorter than this.
#interleaved=auto      If true, forces fastq input to be paired and interleaved.
#ziplevel=2            Set to 1 (lowest) through 9 (max) to change compression level; lower compression is faster.
#
#Output format parameters:
#storename=t           (sn) Store scaffold names (set false to save memory).
##addpairnum=f         Add .1 and .2 to numeric id of read1 and read2.
#storequality=t        (sq) Store quality values for fastq assemblies (set false to save memory).
#uniquenames=t         (un) Ensure all output scaffolds have unique names.  Uses more memory.
#mergenames=f          When a sequence absorbs another, concatenate their headers.
#mergedelimiter=>      Delimiter between merged headers.  Can be a symbol name like greaterthan.
#numbergraphnodes=t    (ngn) Label dot graph nodes with read numbers rather than read names.
#sort=f                Sort output (otherwise it will be random).  Options:
#                         length:  Sort by length
#                         quality: Sort by quality
#                         name:    Sort by name
#                         id:      Sort by input order
#ascending=f           Sort in ascending order.
#ordered=f             Output sequences in input order.  Equivalent to sort=id ascending.
#renameclusters=f      (rnc) Rename contigs to indicate which cluster they are in.
#printlengthinedges=f  (ple) Print the length of contigs in edges.
#
#Processing parameters:
#absorbrc=t            (arc) Absorb reverse-complements as well as normal orientation.
#absorbmatch=t         (am) Absorb exact matches of contigs.
#absorbcontainment=t   (ac) Absorb full containments of contigs.
##absorboverlap=f      (ao) Absorb (merge) non-contained overlaps of contigs (TODO).
#findoverlap=f         (fo) Find overlaps between contigs (containments and non-containments).  Necessary for clustering.
#uniqueonly=f          (uo) If true, all copies of duplicate reads will be discarded, rather than keeping 1.
#rmn=f                 (requirematchingnames) If true, both names and sequence must match.
#usejni=f              (jni) Do alignments in C code, which is faster, if an edit distance is allowed.
#                      This will require compiling the C code; details are in /jni/README.txt.
#
#Subset parameters:
#subsetcount=1         (sstc) Number of subsets used to process the data; higher uses less memory.
#subset=0              (sst) Only process reads whose ((ID%subsetcount)==subset).
#
#Clustering parameters:
#cluster=f             (c) Group overlapping contigs into clusters.
#pto=f                 (preventtransitiveoverlaps) Do not look for new edges between nodes in the same cluster.
#minclustersize=1      (mcs) Do not output clusters smaller than this.
#pbr=f                 (pickbestrepresentative) Only output the single highest-quality read per cluster.
#
#Cluster postprocessing parameters:
#processclusters=f     (pc) Run the cluster processing phase, which performs the selected operations in this category.
#                      For example, pc AND cc must be enabled to perform cc.
#fixmultijoins=t       (fmj) Remove redundant overlaps between the same two contigs.
#removecycles=t        (rc) Remove all cycles so clusters form trees.
#cc=t                  (canonicizeclusters) Flip contigs so clusters have a single orientation.
#fcc=f                 (fixcanoncontradictions) Truncate graph at nodes with canonization disputes.
#foc=f                 (fixoffsetcontradictions) Truncate graph at nodes with offset disputes.
#mst=f                 (maxspanningtree) Remove cyclic edges, leaving only the longest edges that form a tree.
#
#Overlap Detection Parameters
#exact=t               (ex) Only allow exact symbol matches.  When false, an 'N' will match any symbol.
#touppercase=t         (tuc) Convert input bases to upper-case; otherwise, lower-case will not match.
#maxsubs=0             (s) Allow up to this many mismatches (substitutions only, no indels).  May be set higher than maxedits.
#maxedits=0            (e) Allow up to this many edits (subs or indels).  Higher is slower.
#minidentity=100       (mid) Absorb contained sequences with percent identity of at least this (includes indels).
#minlengthpercent=0    (mlp) Smaller contig must be at least this percent of larger contig's length to be absorbed.
#minoverlappercent=0   (mop) Overlap must be at least this percent of smaller contig's length to cluster and merge.
#minoverlap=200        (mo) Overlap must be at least this long to cluster and merge.
#depthratio=0          (dr) When non-zero, overlaps will only be formed between reads with a depth ratio of at most this.
#                      Should be above 1.  Depth is determined by parsing the read names; this information can be added
#                      by running KmerNormalize (khist.sh, bbnorm.sh, or ecc.sh) with the flag 'rename'
#k=31                  Seed length used for finding containments and overlaps.  Anything shorter than k will not be found.
#numaffixmaps=1        (nam) Number of prefixes/suffixes to index per contig. Higher is more sensitive, if edits are allowed.
#hashns=f              Set to true to search for matches using kmers containing Ns.  Can lead to extreme slowdown in some cases.
##ignoreaffix1=f       (ia1) Ignore first affix (for testing).
##storesuffix=f        (ss) Store suffix as well as prefix.  Automatically set to true when doing inexact matches.
#
#Other Parameters
#qtrim=f               Set to qtrim=rl to trim leading and trailing Ns.
#trimq=6               Quality trim level.
#forcetrimleft=-1      (ftl) If positive, trim bases to the left of this position (exclusive, 0-based).
#forcetrimright=-1     (ftr) If positive, trim bases to the right of this position (exclusive, 0-based).
#
#Note on Proteins / Amino Acids
#Dedupe supports amino acid space via the 'amino' flag.  This also changes the default kmer length to 10.
#In amino acid mode, all flags related to canonicity and reverse-complementation are disabled,
#and nam (numaffixmaps) is currently limited to 2 per tip.
#
#Java Parameters:
#-Xmx                  This will set Java's memory usage, overriding autodetection.
#                      -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.
#                    The max is typically 85% of physical memory.
#-eoom                 This flag will cause the process to exit if an out-of-memory exception occurs.  Requires Java 8u92+.
#-da                   Disable assertions.
#
#Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems.
#
#####

##### fastqc #####
#
#            FastQC - A high throughput sequence QC analysis tool
#
#SYNOPSIS
#
#        fastqc seqfile1 seqfile2 .. seqfileN
#
#    fastqc [-o output dir] [--(no)extract] [-f fastq|bam|sam]
#           [-c contaminant file] seqfile1 .. seqfileN
#
#DESCRIPTION
#
#    FastQC reads a set of sequence files and produces from each one a quality
#    control report consisting of a number of different modules, each one of
#    which will help to identify a different potential type of problem in your
#    data.
#
#    If no files to process are specified on the command line then the program
#    will start as an interactive graphical application.  If files are provided
#    on the command line then the program will run with no user interaction
#    required.  In this mode it is suitable for inclusion into a standardised
#    analysis pipeline.
#
#    The options for the program as as follows:
#
#    -h --help       Print this help file and exit
#
#    -v --version    Print the version of the program and exit
#
#    -o --outdir     Create all output files in the specified output directory.
#                    Please note that this directory must exist as the program
#                    will not create it.  If this option is not set then the
#                    output file for each sequence file is created in the same
#                    directory as the sequence file which was processed.
#
#    --casava        Files come from raw casava output. Files in the same sample
#                    group (differing only by the group number) will be analysed
#                    as a set rather than individually. Sequences with the filter
#                    flag set in the header will be excluded from the analysis.
#                    Files must have the same names given to them by casava
#                    (including being gzipped and ending with .gz) otherwise they
#                    won't be grouped together correctly.
#
#    --nano          Files come from nanopore sequences and are in fast5 format. In
#                    this mode you can pass in directories to process and the program
#                    will take in all fast5 files within those directories and produce
#                    a single output file from the sequences found in all files.
#
#    --nofilter      If running with --casava then don't remove read flagged by
#                    casava as poor quality when performing the QC analysis.
#
#    --extract       If set then the zipped output file will be uncompressed in
#                    the same directory after it has been created.  By default
#                    this option will be set if fastqc is run in non-interactive
#                    mode.
#
#    -j --java       Provides the full path to the java binary you want to use to
#                    launch fastqc. If not supplied then java is assumed to be in
#                    your path.
#
#    --noextract     Do not uncompress the output file after creating it.  You
#                    should set this option if you do not wish to uncompress
#                    the output when running in non-interactive mode.
#
#    --nogroup       Disable grouping of bases for reads >50bp. All reports will
#                    show data for every base in the read.  WARNING: Using this
#                    option will cause fastqc to crash and burn if you use it on
#                    really long reads, and your plots may end up a ridiculous size.
#                    You have been warned!
#
#    --min_length    Sets an artificial lower limit on the length of the sequence
#                    to be shown in the report.  As long as you set this to a value
#                    greater or equal to your longest read length then this will be
#                    the sequence length used to create your read groups.  This can
#                    be useful for making directly comaparable statistics from
#                    datasets with somewhat variable read lengths.
#
#    -f --format     Bypasses the normal sequence file format detection and
#                    forces the program to use the specified format.  Valid
#                    formats are bam,sam,bam_mapped,sam_mapped and fastq
#
#    -t --threads    Specifies the number of files which can be processed
#                    simultaneously.  Each thread will be allocated 250MB of
#                    memory so you shouldn't run more threads than your
#                    available memory will cope with, and not more than
#                    6 threads on a 32 bit machine
#
#    -c              Specifies a non-default file which contains the list of
#    --contaminants  contaminants to screen overrepresented sequences against.
#                    The file must contain sets of named contaminants in the
#                    form name[tab]sequence.  Lines prefixed with a hash will
#                    be ignored.
#
#    -a              Specifies a non-default file which contains the list of
#    --adapters      adapter sequences which will be explicity searched against
#                    the library. The file must contain sets of named adapters
#                    in the form name[tab]sequence.  Lines prefixed with a hash
#                    will be ignored.
#
#    -l              Specifies a non-default file which contains a set of criteria
#    --limits        which will be used to determine the warn/error limits for the
#                    various modules.  This file can also be used to selectively
#                    remove some modules from the output all together.  The format
#                    needs to mirror the default limits.txt file found in the
#                    Configuration folder.
#
#   -k --kmers       Specifies the length of Kmer to look for in the Kmer content
#                    module. Specified Kmer length must be between 2 and 10. Default
#                    length is 7 if not specified.
#
#   -q --quiet       Supress all progress messages on stdout and only report errors.
#
#   -d --dir         Selects a directory to be used for temporary files written when
#                    generating report images. Defaults to system temp directory if
#                    not specified.
#
#BUGS
#
#    Any bugs in fastqc should be reported either to simon.andrews@babraham.ac.uk
#    or in www.bioinformatics.babraham.ac.uk/bugzilla/
#
#####










