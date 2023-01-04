#!/bin/bash

##  use mapping strategy for consensus
##  try bbmap with long read params for this

##  requirements
#
#	bbmap
#	rename
#	samtools








##  activate the environment for this downstream analysis
eval "$(conda shell.bash hook)";
conda activate virus-nanopore;

## input files/directories 
DBdir=/mnt/lely_DB/AI-nanopore; 
DB=VI.gisaid.dedupe.fasta;
DBin="$DBdir"/"$DB";
READdir="$PWD"/RAWREADS;
MAPPINGdir="$PWD"/02_mapping;
LOGS="$PWD"/LOGS;
REFdir="$PWD"/references;


## create directory
mkdir -p "$MAPPINGdir";

## counter
count0=1;
countS=$(cat samples.txt | wc -l);

while [ "$count0" -le "$countS" ]; do

	SAMPLE=$(cat samples.txt | awk 'NR=='$count0);

		echo $SAMPLE;
		
#  I/O files
		FileIn="$READdir"/"$SAMPLE".chopped.fastq.gz;
		TEMP="$READdir"/"$SAMPLE".chopped.short.head.fastq.gz;
		FileOut="$MAPPINGdir"/"$SAMPLE".sam;
		SCAFSTATS="$MAPPINGdir"/"$SAMPLE".refstats.tab;
		LOG="$LOGS"/"$SAMPLE".select.ref.mapping.log;


			reformat.sh -Xmx20g in="$FileIn" out="$TEMP" trimreaddescription=t ow;

			mapPacBio.sh -Xmx20g in="$TEMP" out="$FileOut" ref="$DBin" nodisk scafstats="$SCAFSTATS" ambiguous=best qin=33 ow > "$LOG" 2>&1 ; ## only for selecting the proper reference

count0=$((count0+1));
done


##  select the proper reference from the used database


## counter
count1=1;
countT=$(cat samples.txt | wc -l);

while [ "$count1" -le "$countT" ]; do

	sample=$(cat samples.txt | awk 'NR=='"$count1");

		echo $sample;
		
		rm "$MAPPINGdir"/"$sample".ref.name.lst;
		
			for i in PB1 PB2 PA NA HA MP NP NS;do
			
			echo $i;
			
				FileIn="$MAPPINGdir"/"$sample".refstats.tab;
			
					cat "$FileIn" | grep "|$i|" | sort -k9,9 -nr | head -n1 | cut -f1 >> "$MAPPINGdir"/"$sample".ref.name.lst;
					
			done


REFout="$REFdir"/"$sample".selected.refs.fa;


filterbyname.sh in="$DBin" out="$REFout" names="$MAPPINGdir"/"$sample".ref.name.lst include=t fastawrap=10000 ow ; ## select the proper ref sequences for constructing consensus

count1=$((count1+1));
done


##  mapping reads on proper refs and constructing consensus sequence


## counter
count2=1;
countU=$(cat samples.txt | wc -l);

while [ "$count2" -le "$countU" ]; do

	SAMPLE=$(cat samples.txt | awk 'NR=='$count2);

		echo $SAMPLE;
		
#  I/O files
		FileIn="$READdir"/"$SAMPLE".chopped.short.head.fastq.gz;
		FileOut="$MAPPINGdir"/"$SAMPLE".reference.sam;
		FileOutm="$MAPPINGdir"/"$SAMPLE".reference.minimap2.sam;
		FileOut2="$MAPPINGdir"/"$SAMPLE".reference.bam;
		FileOut2m="$MAPPINGdir"/"$SAMPLE".reference.minimap2.bam;		
		REFin="$REFdir"/"$sample".selected.refs.fa;
		CONSoutTemp="$MAPPINGdir"/"$SAMPLE".reference.consensus.temp.fa;
		CONSoutTempm="$MAPPINGdir"/"$SAMPLE".reference.consensus.minimap2.temp.fa;
		CONSout="$MAPPINGdir"/"$SAMPLE".reference.consensus.fa;
		CONSoutm="$MAPPINGdir"/"$SAMPLE".reference.consensus.minimap2.fa;
		LOG="$LOGS"/"$SAMPLE".ref.mapping.log;

# mapping and create consensus
			mapPacBio.sh -Xmx20g in="$FileIn" out="$FileOut" ref="$REFin" nodisk ambiguous=best sam=1.3 qin=33 ow >> "$LOG" 2>&1 ; ## only for selecting the proper reference

			samtools view -bShu "$FileOut" | samtools sort -m 2G -@ 3 - -o "$FileOut2"; 
			samtools index "$FileOut2";		
			samtools consensus -@ 48 -a -q -f fasta --show-ins yes --show-del yes "$FileOut2"  -o "$CONSoutTemp";

cat "$CONSoutTemp" | sed "s/^>/>"$sample"_/g" > "$CONSout";

			
			minimap2 -ax map-ont "$REFin" "$FileIn" > "$FileOutm";

			samtools view -bShu "$FileOutm" | samtools sort -m 2G -@ 3 - -o "$FileOut2m"; 
			samtools index "$FileOut2m";		
			samtools consensus -@ 48 -a -q -f fasta --show-ins yes --show-del yes "$FileOut2m"  -o "$CONSoutTempm";

cat "$CONSoutTempm" | sed "s/^>/>"$sample"_/g" > "$CONSoutm";






count2=$((count2+1));
done











exit 1


##### mapPacBio
#BBMap
#Written by Brian Bushnell, from Dec. 2010 - present
#Last modified September 15, 2022
#
#Description:  Fast and accurate splice-aware read aligner.
#Please read bbmap/docs/guides/BBMapGuide.txt for more information.
#
#To index:     bbmap.sh ref=<reference fasta>
#To map:       bbmap.sh in=<reads> out=<output sam>
#To map without writing an index:
#    bbmap.sh ref=<reference fasta> in=<reads> out=<output sam> nodisk
#
#in=stdin will accept reads from standard in, and out=stdout will write to
#standard out, but file extensions are still needed to specify the format of the
#input and output files e.g. in=stdin.fa.gz will read gzipped fasta from
#standard in; out=stdout.sam.gz will write gzipped sam.
#
#Indexing Parameters (required when building the index):
#nodisk=f                Set to true to build index in memory and write nothing
#                        to disk except output.
#ref=<file>              Specify the reference sequence.  Only do this ONCE,
#                        when building the index (unless using 'nodisk').
#build=1                 If multiple references are indexed in the same directory,
#                        each needs a unique numeric ID (unless using 'nodisk').
#k=13                    Kmer length, range 8-15.  Longer is faster but uses
#                        more memory.  Shorter is more sensitive.
#                        If indexing and mapping are done in two steps, K should
#                        be specified each time.
#path=<.>                Specify the location to write the index, if you don't
#                        want it in the current working directory.
#usemodulo=f             Throw away ~80% of kmers based on remainder modulo a
#                        number (reduces RAM by 50% and sensitivity slightly).
#                        Should be enabled both when building the index AND
#                        when mapping.
#rebuild=f               Force a rebuild of the index (ref= should be set).
#
#Input Parameters:
#build=1                 Designate index to use.  Corresponds to the number
#                        specified when building the index.
#in=<file>               Primary reads input; required parameter.
#in2=<file>              For paired reads in two files.
#interleaved=auto        True forces paired/interleaved input; false forces
#                        single-ended mapping. If not specified, interleaved
#                        status will be autodetected from read names.
#fastareadlen=500        Break up FASTA reads longer than this.  Max is 500 for
#                        BBMap and 6000 for BBMapPacBio.  Only works for FASTA
#                        input (use 'maxlen' for FASTQ input).  The default for
#                        bbmap.sh is 500, and for mapPacBio.sh is 6000.
#unpigz=f                Spawn a pigz (parallel gzip) process for faster
#                        decompression than using Java.
#                        Requires pigz to be installed.
#touppercase=t           (tuc) Convert lowercase letters in reads to upper case
#                        (otherwise they will not match the reference).
#
#Sampling Parameters:
#
#reads=-1                Set to a positive number N to only process the first N
#                        reads (or pairs), then quit.  -1 means use all reads.
#samplerate=1            Set to a number from 0 to 1 to randomly select that
#                        fraction of reads for mapping. 1 uses all reads.
#skipreads=0             Set to a number N to skip the first N reads (or pairs),
#                        then map the rest.
#
#Mapping Parameters:
#fast=f                  This flag is a macro which sets other paramters to run
#                        faster, at reduced sensitivity.  Bad for RNA-seq.
#slow=f                  This flag is a macro which sets other paramters to run
#                        slower, at greater sensitivity.  'vslow' is even slower.
#maxindel=16000          Don't look for indels longer than this. Lower is faster.
#                        Set to >=100k for RNAseq with long introns like mammals.
#strictmaxindel=f        When enabled, do not allow indels longer than 'maxindel'.
#                        By default these are not sought, but may be found anyway.
#tipsearch=100           Look this far for read-end deletions with anchors
#                        shorter than K, using brute force.
#minid=0.76              Approximate minimum alignment identity to look for.
#                        Higher is faster and less sensitive.
#minhits=1               Minimum number of seed hits required for candidate sites.
#                        Higher is faster.
#local=f                 Set to true to use local, rather than global, alignments.
#                        This will soft-clip ugly ends of poor alignments.
#perfectmode=f           Allow only perfect mappings when set to true (very fast).
#semiperfectmode=f       Allow only perfect and semiperfect (perfect except for
#                        N's in the reference) mappings.
#threads=auto            (t) Set to number of threads desired.  By default, uses
#                        all cores available.
#ambiguous=best          (ambig) Set behavior on ambiguously-mapped reads (with
#                        multiple top-scoring mapping locations).
#                            best    (use the first best site)
#                            toss    (consider unmapped)
#                            random  (select one top-scoring site randomly)
#                            all     (retain all top-scoring sites)
#samestrandpairs=f       (ssp) Specify whether paired reads should map to the
#                        same strand or opposite strands.
#requirecorrectstrand=t  (rcs) Forbid pairing of reads without correct strand
#                        orientation.  Set to false for long-mate-pair libraries.
#killbadpairs=f          (kbp) If a read pair is mapped with an inappropriate
#                        insert size or orientation, the read with the lower
#                        mapping quality is marked unmapped.
#pairedonly=f            (po) Treat unpaired reads as unmapped.  Thus they will
#                        be sent to 'outu' but not 'outm'.
#rcomp=f                 Reverse complement both reads prior to mapping (for LMP
#                        outward-facing libraries).
#rcompmate=f             Reverse complement read2 prior to mapping.
#pairlen=32000           Set max allowed distance between paired reads.
#                        (insert size)=(pairlen)+(read1 length)+(read2 length)
#rescuedist=1200         Don't try to rescue paired reads if avg. insert size
#                        greater than this.  Lower is faster.
#rescuemismatches=32     Maximum mismatches allowed in a rescued read.  Lower
#                        is faster.
#averagepairdist=100     (apd) Initial average distance between paired reads.
#                        Varies dynamically; does not need to be specified.
#deterministic=f         Run in deterministic mode.  In this case it is good
#                        to set averagepairdist.  BBMap is deterministic
#                        without this flag if using single-ended reads,
#                        or run singlethreaded.
#bandwidthratio=0        (bwr) If above zero, restrict alignment band to this
#                        fraction of read length.  Faster but less accurate.
#bandwidth=0             (bw) Set the bandwidth directly.
#                        fraction of read length.  Faster but less accurate.
#usejni=f                (jni) Do alignments faster, in C code.  Requires
#                        compiling the C code; details are in /jni/README.txt.
#maxsites2=800           Don't analyze (or print) more than this many alignments
#                        per read.
#ignorefrequentkmers=t   (ifk) Discard low-information kmers that occur often.
#excludefraction=0.03    (ef) Fraction of kmers to ignore.  For example, 0.03
#                        will ignore the most common 3% of kmers.
#greedy=t                Use a greedy algorithm to discard the least-useful
#                        kmers on a per-read basis.
#kfilter=0               If positive, potential mapping sites must have at
#                        least this many consecutive exact matches.
#
#
#Quality and Trimming Parameters:
#qin=auto                Set to 33 or 64 to specify input quality value ASCII
#                        offset. 33 is Sanger, 64 is old Solexa.
#qout=auto               Set to 33 or 64 to specify output quality value ASCII
#                        offset (only if output format is fastq).
#qtrim=f                 Quality-trim ends before mapping.  Options are:
#                        'f' (false), 'l' (left), 'r' (right), and 'lr' (both).
#untrim=f                Undo trimming after mapping.  Untrimmed bases will be
#                        soft-clipped in cigar strings.
#trimq=6                 Trim regions with average quality below this
#                        (phred algorithm).
#mintrimlength=60        (mintl) Don't trim reads to be shorter than this.
#fakefastaquality=-1     (ffq) Set to a positive number 1-50 to generate fake
#                        quality strings for fasta input reads.
#ignorebadquality=f      (ibq) Keep going, rather than crashing, if a read has
#                        out-of-range quality values.
#usequality=t            Use quality scores when determining which read kmers
#                        to use as seeds.
#minaveragequality=0     (maq) Do not map reads with average quality below this.
#maqb=0                  If positive, calculate maq from this many initial bases.
#
#Output Parameters:
#out=<file>              Write all reads to this file.
#outu=<file>             Write only unmapped reads to this file.  Does not
#                        include unmapped paired reads with a mapped mate.
#outm=<file>             Write only mapped reads to this file.  Includes
#                        unmapped paired reads with a mapped mate.
#mappedonly=f            If true, treats 'out' like 'outm'.
#bamscript=<file>        (bs) Write a shell script to <file> that will turn
#                        the sam output into a sorted, indexed bam file.
#ordered=f               Set to true to output reads in same order as input.
#                        Slower and uses more memory.
#overwrite=f             (ow) Allow process to overwrite existing files.
#secondary=f             Print secondary alignments.
#sssr=0.95               (secondarysitescoreratio) Print only secondary alignments
#                        with score of at least this fraction of primary.
#ssao=f                  (secondarysiteasambiguousonly) Only print secondary
#                        alignments for ambiguously-mapped reads.
#maxsites=5              Maximum number of total alignments to print per read.
#                        Only relevant when secondary=t.
#quickmatch=f            Generate cigar strings more quickly.
#trimreaddescriptions=f  (trd) Truncate read and ref names at the first whitespace,
#                        assuming that the remainder is a comment or description.
#ziplevel=2              (zl) Compression level for zip or gzip output.
#pigz=f                  Spawn a pigz (parallel gzip) process for faster
#                        compression than Java.  Requires pigz to be installed.
#machineout=f            Set to true to output statistics in machine-friendly
#                        'key=value' format.
#printunmappedcount=f    Print the total number of unmapped reads and bases.
#                        If input is paired, the number will be of pairs
#                        for which both reads are unmapped.
#showprogress=0          If positive, print a '.' every X reads.
#showprogress2=0         If positive, print the number of seconds since the
#                        last progress update (instead of a '.').
#renamebyinsert=f        Renames reads based on their mapped insert size.
#
#Bloom-Filtering Parameters (bloomfilter.sh is the standalone version).
#bloom=f                 Use a Bloom filter to ignore reads not sharing kmers
#                        with the reference.  This uses more memory, but speeds
#                        mapping when most reads don't match the reference.
#bloomhashes=2           Number of hash functions.
#bloomminhits=3          Number of consecutive hits to be considered matched.
#bloomk=31               Bloom filter kmer length.
#bloomserial=t           Use the serialized Bloom filter for greater loading
#                        speed, if available.  If not, generate and write one.
#
#Post-Filtering Parameters:
#idfilter=0              Independant of minid; sets exact minimum identity
#                        allowed for alignments to be printed.  Range 0 to 1.
#subfilter=-1            Ban alignments with more than this many substitutions.
#insfilter=-1            Ban alignments with more than this many insertions.
#delfilter=-1            Ban alignments with more than this many deletions.
#indelfilter=-1          Ban alignments with more than this many indels.
#editfilter=-1           Ban alignments with more than this many edits.
#inslenfilter=-1         Ban alignments with an insertion longer than this.
#dellenfilter=-1         Ban alignments with a deletion longer than this.
#nfilter=-1              Ban alignments with more than this many ns.  This
#                        includes nocall, noref, and off scaffold ends.
#
#Sam flags and settings:
#noheader=f              Disable generation of header lines.
#sam=1.4                 Set to 1.4 to write Sam version 1.4 cigar strings,
#                        with = and X, or 1.3 to use M.
#saa=t                   (secondaryalignmentasterisks) Use asterisks instead of
#                        bases for sam secondary alignments.
#cigar=t                 Set to 'f' to skip generation of cigar strings (faster).
#keepnames=f             Keep original names of paired reads, rather than
#                        ensuring both reads have the same name.
#intronlen=999999999     Set to a lower number like 10 to change 'D' to 'N' in
#                        cigar strings for deletions of at least that length.
#rgid=                   Set readgroup ID.  All other readgroup fields
#                        can be set similarly, with the flag rgXX=
#                        If you set a readgroup flag to the word 'filename',
#                        e.g. rgid=filename, the input file name will be used.
#mdtag=f                 Write MD tags.
#nhtag=f                 Write NH tags.
#xmtag=f                 Write XM tags (may only work correctly with ambig=all).
#amtag=f                 Write AM tags.
#nmtag=f                 Write NM tags.
#xstag=f                 Set to 'xs=fs', 'xs=ss', or 'xs=us' to write XS tags
#                        for RNAseq using firststrand, secondstrand, or
#                        unstranded libraries.  Needed by Cufflinks.
#                        JGI mainly uses 'firststrand'.
#stoptag=f               Write a tag indicating read stop location, prefixed by YS:i:
#lengthtag=f             Write a tag indicating (query,ref) alignment lengths,
#                        prefixed by YL:Z:
#idtag=f                 Write a tag indicating percent identity, prefixed by YI:f:
#inserttag=f             Write a tag indicating insert size, prefixed by X8:Z:
#scoretag=f              Write a tag indicating BBMap's raw score, prefixed by YR:i:
#timetag=f               Write a tag indicating this read's mapping time, prefixed by X0:i:
#boundstag=f             Write a tag indicating whether either read in the pair
#                        goes off the end of the reference, prefixed by XB:Z:
#notags=f                Turn off all optional tags.
#
#Histogram and statistics output parameters:
#scafstats=<file>        Statistics on how many reads mapped to which scaffold.
#refstats=<file>         Statistics on how many reads mapped to which reference
#                        file; only for BBSplit.
#sortscafs=t             Sort scaffolds or references by read count.
#bhist=<file>            Base composition histogram by position.
#qhist=<file>            Quality histogram by position.
#aqhist=<file>           Histogram of average read quality.
#bqhist=<file>           Quality histogram designed for box plots.
#lhist=<file>            Read length histogram.
#ihist=<file>            Write histogram of insert sizes (for paired reads).
#ehist=<file>            Errors-per-read histogram.
#qahist=<file>           Quality accuracy histogram of error rates versus
#                        quality score.
#indelhist=<file>        Indel length histogram.
#mhist=<file>            Histogram of match, sub, del, and ins rates by
#                        read location.
#gchist=<file>           Read GC content histogram.
#gcbins=100              Number gchist bins.  Set to 'auto' to use read length.
#gcpairs=t               Use average GC of paired reads.
#idhist=<file>           Histogram of read count versus percent identity.
#idbins=100              Number idhist bins.  Set to 'auto' to use read length.
#statsfile=stderr        Mapping statistics are printed here.
#
#Coverage output parameters (these may reduce speed and use more RAM):
#covstats=<file>         Per-scaffold coverage info.
#rpkm=<file>             Per-scaffold RPKM/FPKM counts.
#covhist=<file>          Histogram of # occurrences of each depth level.
#basecov=<file>          Coverage per base location.
#bincov=<file>           Print binned coverage per location (one line per X bases).
#covbinsize=1000         Set the binsize for binned coverage output.
#nzo=t                   Only print scaffolds with nonzero coverage.
#twocolumn=f             Change to true to print only ID and Avg_fold instead of
#                        all 6 columns to the 'out=' file.
#32bit=f                 Set to true if you need per-base coverage over 64k.
#strandedcov=f           Track coverage for plus and minus strand independently.
#startcov=f              Only track start positions of reads.
#secondarycov=t          Include coverage of secondary alignments.
#physcov=f               Calculate physical coverage for paired reads.
#                        This includes the unsequenced bases.
#delcoverage=t           (delcov) Count bases covered by deletions as covered.
#                        True is faster than false.
#covk=0                  If positive, calculate kmer coverage statistics.
#
#Java Parameters:
#-Xmx                    This will set Java's memory usage,
#                        overriding autodetection.
#                        -Xmx20g will specify 20 gigs of RAM, and -Xmx800m
#                        will specify 800 megs.  The max is typically 85% of
#                        physical memory.  The human genome requires around 24g,
#                        or 12g with the 'usemodulo' flag.  The index uses
#                        roughly 6 bytes per reference base.
#-eoom                   This flag will cause the process to exit if an
#                        out-of-memory exception occurs.  Requires Java 8u92+.
#-da                     Disable assertions.
#
#Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter
#any problems, or post at: http://seqanswers.com/forums/showthread.php?t=41057
#
#####

##### filterbyname
#
#Written by Brian Bushnell
#Last modified September 1, 2016
#
#Description:  Filters reads by name.
#
#Usage:  filterbyname.sh in=<file> in2=<file2> out=<outfile> out2=<outfile2> names=<string,string,string> include=<t/f>
#
#in2 and out2 are for paired reads and are optional.
#If input is paired and there is only one output file, it will be written interleaved.
#Important!  Leading > and @ symbols are NOT part of sequence names;  they are part of
#the fasta, fastq, and sam specifications.  Therefore, this is correct:
#names=e.coli_K12
#And these are incorrect:
#names=>e.coli_K12
#names=@e.coli_K12
#
#Parameters:
#include=f       Set to 'true' to include the filtered names rather than excluding them.
#substring=f     Allow one name to be a substring of the other, rather than a full match.
#                   f: No substring matching.
#                   t: Bidirectional substring matching.
#                   header: Allow input read headers to be substrings of names in list.
#                   name: Allow names in list to be substrings of input read headers.
#prefix=f        Allow names to match read header prefixes.
#case=t          (casesensitive) Match case also.
#ow=t            (overwrite) Overwrites files that already exist.
#app=f           (append) Append to files that already exist.
#zl=4            (ziplevel) Set compression level, 1 (low) to 9 (max).
#int=f           (interleaved) Determines whether INPUT file is considered interleaved.
#names=          A list of strings or files.  The files can have one name per line, or
#                be a standard read file (fasta, fastq, or sam).
#minlen=0        Do not output reads shorter than this.
#ths=f           (truncateheadersymbol) Ignore a leading @ or > symbol in the names file.
#tws=f           (truncatewhitespace) Ignore leading or trailing whitespace in the names file.
#truncate=f      Set both ths and tws at the same time.
#
#Positional parameters:
#These optionally allow you to output only a portion of a sequence.  Zero-based, inclusive.
#Intended for a single sequence and include=t mode.
#from=-1         Only print bases starting at this position.
#to=-1           Only print bases up to this position.
#range=          Set from and to with a single flag.
#
#
#Java Parameters:
#-Xmx            This will set Java's memory usage, overriding autodetection.
#                -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.
#                    The max is typically 85% of physical memory.
#-eoom           This flag will cause the process to exit if an out-of-memory
#                exception occurs.  Requires Java 8u92+.
#-da             Disable assertions.
#
#To read from stdin, set 'in=stdin'.  The format should be specified with an extension, like 'in=stdin.fq.gz'
#To write to stdout, set 'out=stdout'.  The format should be specified with an extension, like 'out=stdout.fasta'
#
#Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems.
#
#####
