#!/usr/bin/perl -w

use strict;
use File::Copy;
use Getopt::Std;
use POSIX;
use GTF;
use FAlite;
use BPlite;
use EGlite;
use Sim4Parser;


use vars qw($opt_p $opt_r $opt_h $opt_c $opt_s $opt_t $opt_n $opt_o $opt_m); 
getopts('pr:hcs:tznom');

my $usage = "usage: $0 [opts] <parameterfile> <fragment>
OR -m <parameterfile> <gtf file>

This program uses two methods of finding pseudogenes in input GTF files, the intron alignment method and the conserved synteny method (described by Van Baren and Brent, manuscript in preparation). The program can run on 1 Mb fragments as used by Twinscan, or on whole chromosome GTFs using the m option.
In its default state, it takes in a fragment directory, runs the pseudogene finding pipeline on the gtf in that fragment (named N.seq.gtf), masks the sequence in the genomic sequence found in the fragment (named N.seq.masked and reruns Twinscan. This process is repeated until no more pseudogenes are found. In whole chromosome mode, it takes the chromosome as argument of the m option, and a full path to a chromosome GTF file. For GTF format, see website. In this mode, it outputs a list of genes that contain pseudogene exons, as well as a GTF file annotating those pseudogene regions.
  -h: Display this help message
  -t: do not delete files in the temporary directory (for editing program only)
  -r <nr> : number of rounds to go through. Default is 6 (including intron method).
  -o: Do not run the intron alignment method (after the protein method)
  -p: Do not run the protein method (go to intron alignment method directly)
  -c: After first gene prediction run, output and exit;
  -n: Stop after masking, do not run Twinscan
  -m: Find pseudogenes then stop. In this case, do input the gtf file named chrN.(*)gtf and its FULL PATH!

";

# The program takes a directory as input, check if it exists
our($paramfile, $homedir) = @ARGV;
die $usage unless ($homedir);
if($opt_h){
    print $usage;
    exit(0);
}
die "cannot find $homedir\n" unless ( -e $homedir );

# get all parameters from parameter file
our($BLAST, $BLASTP, $BLASTN, $TBLASTN, $SIM4, $EST2GENO, $parent_blastfile, $parent_geno_seqs, $parent_clone_seqs, 
	$target_seqdir, $informant_seqdir, $prot_db, $prot_locations_file, $protein_sequences, $prot_loc_inputdir, 
	$synteny_dir, $localdir, $opt_m, $species, $alignment_program, $align_ext, $align_path,
	$informant_blastdb, $nscan_params, $iscan_params, $TWINSCAN_COMMAND, $run_gene_pred);

my@errors = format_params($paramfile);

if(exists ($errors[0])){
	foreach(@errors){
		print "$_\n";
	}
	print "quitting...\n";
	exit;
}

# create a comment array (output in pgene.reason on exit)
# must be accessible all through program
our@comments = ();

# the coords hash will hold the protein locations  
our%coords;

# Program only runs maximum 9 rounds of iterative masking
# (gene prediction is less often)
# If this is exceeded, something is probably wrong with
# the code and/or the input
our$round = "1";
my$round_cutoff = $round + 9;

# create a temporary outputdir for all temporary files
# seed rand
srand(time() ^ ($$ + ($$ << 15)));
my$nr = rand;
$nr =~ s/^0\.//;

my$tmpdir="/$localdir/tmp.$nr";
graceful_die ("Unique dir $tmpdir exists") if ( -d "$tmpdir");
mkdir "$tmpdir" || graceful_die ("cannot create temporary directory $tmpdir: $!");

# Create a temporary directory for files that are going to be copied back to the working directory
our$returndir = "$tmpdir/return";
mkdir "$returndir" || graceful_die ("cannot create temporary directory $returndir: $!");

# Get the chromosome and fragment number when running in default mode
# When running whole chromosomes, get info from input arguments
# also create clone_frag hash for filling when $opt_m is run (whole chromosomes). This keeps
# track of the location of genes to compare with predicted genes
our($chr, $frag);
our%clone_frag=();

if($opt_m){
  ($chr)=$homedir=~/chr([^\/]*?)\.[^\/]*gtf$/;
  die "for m option, input should be a file named chrN(.*).gtf!\n" unless $chr;
  $frag="$chr";
}else{
  ($chr, $frag) = $homedir =~ /chr_split\/(.*?)\/fragment_(\d+)/;
  graceful_die ("ERROR inputdir $homedir does not seem to be a fragment") unless (defined($frag) && defined($chr));
  graceful_die ("not an existing directory") unless ( -e "$homedir");
  graceful_die ("$homedir does not look like a full path") unless ($homedir =~ /^\//);
}

# only work on fragments that contain a gtf file (when running gene prediction)
our$gff_file = "$homedir\/$frag.seq.gtf";

# if the gtf file is very big (as expected in opt_m files), split per ~5000 genes and run for
# each sub gtf separately. This means spawning children and collecting the output.
if($opt_m){
  die "$homedir should be full path!" unless($homedir=~ /^\//);
  copy ("$homedir", "$returndir/") || graceful_die ("cannot copy $homedir to $returndir: $!");
  $homedir=~ s/([\w\d\.]+?)$//;
  $gff_file=$&;
  my@outfiles=split_gtf_for_optm("$returndir/$gff_file");

# if the first entry in @outfiles is the gff file, it has not been split. If it has been split, 
# create child process for every file
  if($outfiles[0] !~ /$gff_file$/){
# remove previously created tmpfile
    system "rm -r $tmpdir";
    foreach my$childfile(@outfiles){
      if($opt_p){
        system "$0 -p -m $paramfile $childfile";
      }elsif($opt_o){
        system "$0 -o -m $paramfile $childfile";
      }else{
        system "$0 -m $paramfile $childfile";
      }

# the system will copy the files to the 'homedir' -the place from which the run was started
# remove temporary directories and concatenate outputfiles
       (my$childdir=$childfile)=~ s/\/[\w\d\.]+$//;
       if( -e "$childdir/$frag.ERROR"){
         system "cat $childdir/$frag.ERROR >> $homedir/$frag.ERROR";
       }else{
         system "cat $childdir/*.masked.txt >> $homedir/$frag.masked.txt";
         system "cat $childdir/$frag.masked.gtf >> $homedir/$frag.masked.gtf";
         system "cat $childdir/$chr.pgene.reason >> $homedir/$frag.pgene.reason";
         system "rm -r $childdir";
       }
    }
    exit;
  }
  $gff_file="$homedir/$gff_file";
}

unless ( -e "$gff_file"){
  graceful_exit("no gff file $gff_file in fragment");
}


if($opt_o){
  push(@comments, "run without intron alignment method\n");
}

if($opt_t){
  push(@comments, "found do-not-delete flag\n");
}

# $fragmentseq needed for getting gene prediction sequences. Copy to $returndir
my$fragmentseq;
if($opt_m){
  $fragmentseq="$target_seqdir/chr$chr.masked.fa";
  graceful_die ("cannot find $fragmentseq") unless ( -e $fragmentseq );
  system "ln -s $fragmentseq $tmpdir/$frag.seq.masked";
}else{
  $fragmentseq = "$homedir\/$frag.seq.masked";
  unless ( -e $fragmentseq ){
    $fragmentseq = "$homedir\/$frag.seq";
  }
  copy("$fragmentseq", "$tmpdir/$frag.seq.masked") || 
    graceful_die ("cannot make a local copy of $frag.seq.masked in $tmpdir!: $!");
    $fragmentseq = "$tmpdir/$frag.seq.masked";
}

# only work on fragments haven't already been run
if ( -e "$homedir/$chr.pgene.reason" && !$opt_m){
  graceful_die ("Fragment chr$chr $frag was already done before");
}

# masterseq is the sequence used for pseudogene masking and running Twinscan
# as soon as something is masked, this sequence will change to N.seq.pmasked.masked
my$master_seq = "$tmpdir/$frag.seq.masked";

if( -e  "$homedir/$frag.seq.pmasked.masked" ){
  copy("$homedir/$frag.seq.pmasked.masked", "$returndir/$frag.seq.pmasked.masked");
  $master_seq = "$returndir/$frag.seq.pmasked.masked";
}


# only work on fragments that contain a gff file
unless ( -e "$gff_file"){
  graceful_die("no gff file $gff_file in fragment");
}

# check gff file for entries
my$worthwile;
open(GFF, "$gff_file") || die "cannot open $gff_file: $!\n";
while(<GFF>){
  next if ($_ =~ /^\#|^$/);
  $worthwile="1";
}
graceful_exit ("$gff_file is empty") unless $worthwile;

# our gff file needed in check_introns subroutine. Copy to $returndir and redefine variable
copy ("$gff_file", "$returndir/$frag.seq.gtf") || graceful_die ("cannot copy $gff_file to $returndir: $!");
$gff_file = "$returndir/$frag.seq.gtf";

################ PROGRAM START ###############################

# start the program: run iterative rounds of synteny method and intron alignment method
# (unless one of the methods was excluded by the -p or -o flag)
# Exit after both method if $opt_m is set

my$nochange_in_synteny_method = "0";
my$nochange_in_intron_method = "0";

while(1){
# run synteny method
  unless($opt_p){
    $nochange_in_synteny_method = "0";
    my$current_round = $round;
    my$reason_to_quit = run_synteny_method();
    push(@comments, "protein method: $reason_to_quit\n");
    if($current_round == $round){
      $nochange_in_synteny_method = "1";
    }
    graceful_exit ("Not running intron alignment method: option o set\n$reason_to_quit\n") if ($opt_o);
    $round++;
  }
  push(@comments, "continuing with intron method\n");

# remove everyting created in the temporary directory
  unless($opt_t){
    unlink glob("$tmpdir\/$chr.$frag.*");
  }

# run intron alignment method
  unless($opt_o){
    $nochange_in_intron_method = "0";
    my$current_round = $round;
    my$reason_to_quit = run_intron_method();
    push(@comments, "intron_alignment method: $reason_to_quit\n");
    if($current_round = $round){
      $nochange_in_intron_method = "1";
    }
    graceful_exit ("Not running synteny method: option p set\n$reason_to_quit\n") if ($opt_p);
  }
  if($opt_m){
    graceful_exit("option m set, quitting after intron alignment method\n");
  }
  if($nochange_in_intron_method && $nochange_in_synteny_method){
    graceful_exit ("no masking in either method in round $round\n");
  }
  $round++;
  push(@comments, "continuing with synteny method\n");

# remove everyting created in the temporary directory
  unless($opt_t){
    unlink glob("$tmpdir\/$chr.$frag.*");
  }
}



exit;




### SUBROUTINES ###


# run_synteny_method runs blast against the protein db and compares
# the hits and the query with their genomic locations. If the genomic locations
# do not correspond, the syntenic region in the informant is used to look for
# homology. If no homology is found, the hit is a pseudogene and the corresponding
# region is masked out.
# if $opt_m is set, the subroutine exits before running Twinscan
# called from main program

sub run_synteny_method{

  while (1){

    push(@comments, "starting round $round using synteny method\n");

# get the protein sequences for the predictions  
    my$twinscan_prot = "$tmpdir/$chr.$frag.prot";
    my$protein_seqs = `/bio/bin/gtf2fa.pl -noGTF -frame -ptx -cds -split $gff_file $fragmentseq`;
    if($opt_o){
      graceful_exit ("no protein sequences after translation in round $round") unless $protein_seqs;
    }else{
      return ("no protein sequences after translation in round $round") unless $protein_seqs;
    }

    open (OUTF, ">$twinscan_prot") || graceful_die ("cannot open $twinscan_prot: $!");
    print OUTF $protein_seqs;
    close OUTF;
    remove_small($twinscan_prot);

# Run blast
    my$blastout = "$tmpdir\/$chr.$frag.prot.blastout";
    graceful_die ("cannot access $prot_db") unless ( -e $prot_db );
# The Q parameter is set high to avoid gappy extensions of hits
    system "$BLASTP $prot_db $twinscan_prot -cpus=1 -warnings -notes -E=0.01 -Q=50 > $blastout";

    push(@comments,  "finished first blast\n");

# Get pseudogenes. Output is 
# array with human/informant orthologous regions for the hits
# array with the ID of the twinscan gene with a hit
# array wit the IDs of all the proteins found as parents
    my $synteny_file = "$synteny_dir/synteny_net_chr$chr\.txt";
    graceful_exit("ERROR: $synteny_file not found") unless (-e $synteny_file);
# syntenic region is a hash with all the regions for a given gene
    my ($all_synteny_hits, $hit_ids, $protein_hits, $syntenic_region) 
      = get_pseudogenes_from_protein_blastreport ($blastout, $prot_locations_file, $synteny_file);

# remove the blastoutput file: it can be very big
    unless($opt_t){
      unlink "$blastout";
    }

    if($opt_o){
      graceful_exit("no putative pseudogenes after first blast in round $round") unless ($$all_synteny_hits[0]);
    }else{
      return("no putative pseudogenes after first blast in round $round") unless ($$all_synteny_hits[0]);
    }

# collect the pseudogenes for blast vs informant genome (from the previously created protein file)
    my$informant_blastinput = "/$tmpdir/$chr.$frag.informant_blastinput";
    create_blastinput($twinscan_prot, $informant_blastinput, @$hit_ids);
    
# create blast db from informant synteny sequences (using the array containing the orthologous regions)
    my@informant_coords = get_informant_coords($all_synteny_hits);
    my$informantblastfile = "$tmpdir\/$chr.$frag.informant";
    create_informant_blastdb($informantblastfile, @informant_coords);

# run tblastn vs informant genome fragments
    my$blastoutput="$tmpdir\/$chr.$frag.informant.blastout";
    system "$TBLASTN $informantblastfile $informant_blastinput -warnings -notes -E=0.01 > $blastoutput";

    push (@comments, "ran blast_2\n");

# parse blastoutput. Goodhits contains all informant hits
    my@goodhits = get_good_informant_hits($blastoutput, $syntenic_region);
# save the hits (this file is later put in $homedir as info)
    my$rs_hitfile = "$returndir/$frag.$round.RS_hits";
    open(OUT, ">$rs_hitfile") || graceful_die ("cannot open $rs_hitfile: $!");
    foreach(@$protein_hits){
      print OUT "$_\n";
    }
    print OUT "\n";
    foreach(@goodhits){
      print OUT "$_\n";
    }
    close OUT;


# remove informant hits from pseudogene hits and test for possible nonprocessed pseudogene or paralogs
    my($real_pseudogenes, $info) = compare_blastoutputs($protein_hits, @goodhits);

# mask until oblivion
# this keeps masking the fragment until no more hits are found. The alignment programs have problems
# aligning ests with insertions to genomic sequences. Masking out the first region and realigning
# will find the second region etc. 
    my$mask_flag = "0";
    if($$real_pseudogenes[0]){
      $mask_flag = mask_until_oblivion($real_pseudogenes, "protein");
    }

# create outputfile
    my$pseudogenes_file = "$returndir/$frag.$round.pseudo";
    open (OUTF, ">$pseudogenes_file") || graceful_die ("cannot open outputfile $pseudogenes_file: $!");
    foreach(@$real_pseudogenes){
      print OUTF "$_\n";
    }
    print OUTF "\ninfo:\n";
    foreach(@$info){
      print OUTF "$_\n";
    }
    close OUTF;

    if($opt_m){
      $round++;
      return("Option m set");
    }

    if($opt_o){
      graceful_exit("no pseudogene hits in round $round") unless ($$real_pseudogenes[0]);
      graceful_exit("nothing was masked in round $round") unless ($mask_flag);
    }else{
      return ("no pseudogene hits in round $round") unless ($$real_pseudogenes[0]);
      return ("nothing was masked in round $round") unless ($mask_flag);
    }


    if($opt_n){
      graceful_exit("Option n set: exiting before Twinscan run");
     }

# with the newly masked file, run Twinscan (or nscan)
    run_twinscan();

# add one to counter, and rerun
    $round++;
    if ($round > $round_cutoff ){
      graceful_exit("ERROR: round number too high!");
    }
  }
}

                                                                                

# subtract array removes hits in second array from (referenced) first array

sub subtract_array{
  my($total) = shift(@_);
  my(@subset) = (@_);
  my%hs = ();
  my@rest = ();

#  foreach my$line(@$total){
    foreach my$entry(@subset){
        $hs{$entry}="1";
  }
  foreach my$line(@$total){
    unless(exists$hs{$line}){
      push(@rest, $line);
    }
  }
  return \@rest;
}

# run_twinscan runs nscan, iscan or twinscan, depending on the command in the parameterfile
# called from run_intron_method and run_synteny_method

sub run_twinscan{
	if($run_gene_pred eq "twinscan" || $run_gene_pred eq "iscan"){
		push(@comments, "round $round $TWINSCAN_COMMAND -d $returndir $master_seq\n");
		system "$TWINSCAN_COMMAND -d $returndir $master_seq\n";

# keep only the latest conseq file
		if (-e "$returndir/$frag.conseq"){
			unlink "$returndir/$frag.conseq" || push(@comments, "cannot unlink conseq!");
		}
		move("$master_seq.conseq", "$returndir/$frag.conseq") || push(@comments, "cannot rename conseq!");
# keep the old gff file and put new gff file in place of old one
		move("$gff_file","$returndir\/$frag.$round.gtf");
		move("$master_seq.gtf","$gff_file");

	}elsif($run_gene_pred eq "nscan"){
# copy parameterfile to local directory if it doesn't already exist
		unless ( -e "$tmpdir/$frag.$align_ext"){
			system "cp $align_path/$chr/fragment_$frag/$frag.$align_ext $tmpdir/$frag.$align_ext";
		}
		push(@comments, "round $round $TWINSCAN_COMMAND $master_seq $tmpdir/$frag.$align_ext\n");
		system "$TWINSCAN_COMMAND $master_seq $tmpdir/$frag.$align_ext";
# keep the old gff file and put new gff file in place of old one
		move("$gff_file","$returndir\/$frag.$round.gtf");
		copy("$master_seq.gtf","$gff_file") || die "cannot copy $master_seq.gtf to $gff_file: $!\n";
		unless ( -s $gff_file ){
			graceful_die ("$gff_file is empty, nscan failed");
		}
	}


# remove unnecessary Twinscan files
# remove everyting created in the temporary directory
	unless($opt_t){
		unlink glob("$master_seq.*");
		unlink glob("$tmpdir\/$chr.$frag.*");
	}

# check gff file for entries
	my$worthwile;
	open(GFF, "$gff_file") || die "cannot open $gff_file: $!\n";
	while(<GFF>){
		next if ($_ =~ /^\#|^$/);
		$worthwile="1";
	}
	graceful_exit ("$gff_file is empty") unless $worthwile;
}
                                                                         

# mask_until_oblivion masks the genomic sequence for pseudogenes. After masking
# it retrieves the gene sequence again from the masked file and reruns the
# pseudogene finding process, until no more hits are found.
# called from run_intron_method and run_synteny_method

sub mask_until_oblivion{

  my($real_pseudogenes, $parent) = (@_);
  my$something_was_masked = "0";

# get the protein nt sequences. First, retrieve the chromosome
# the coords hash created in the get_pseudogenes subroutine contains the
# chromosome location of the protein. Extract this chromosome and get the protein sequence
  my$pgene_id_list = "$tmpdir\/$chr.$frag.pgeneid.txt";
  my$estgen_input = "$tmpdir\/$chr.$frag.estgen";
  my@pgid=();

# list for retrieving genomic sequence files
  open (PGID, ">$pgene_id_list") || graceful_die ("cannot open $pgene_id_list: $!");       

# list used for finding combinations of Twinscan genes and proteins
  open (RSLIKE, ">$estgen_input");      
  my@estgen_input= ();

  if($parent eq "protein"){
    foreach my$line(@$real_pseudogenes){
      my($gene, $prot) = $line =~ /^(.*?).exon.*\s.*\s([\w_\.]*)$/;
      $gene=~ s/ multi//;
      my$outfile;

      push(@pgid, $gene);
      print PGID "$gene\n";
      print RSLIKE "$gene\t$prot\n";
      push(@estgen_input, "$gene\t$prot");

      unless($prot =~ /^chr/){
# first see if the protein nt sequence already exists
# the -s flag returns the size of the file
        $outfile = "$protein_sequences/$prot.fa";
        next if (-s "$outfile");
      }

      unless ($coords{$prot}){
        push(@comments, "WARNING cannot find $prot in the protein locations\n");
        next;
      }
      my@protein_positions = split("NEXT", $coords{$prot});
      shift(@protein_positions);

# there may be several chromosome hits per protein, take only the first
# (the sequence will be near identical)
      my($name, $protchr, $rest) = split ("\t", $protein_positions[0]);

# get protein nt sequences, except when using the bootstrapping method 
# (because then we already have them) 
      unless($prot =~ /^chr/){
        my$found = get_protein_sequence($prot, $protchr, $outfile, $prot_loc_inputdir);
        unless($found){
          push(@comments, "WARNING, did not find $prot in subroutine!\n");
        }
      }
    }
  }elsif($parent eq "refseq"){
    foreach my$line(@$real_pseudogenes){
      chomp $line;
      my($gene_id) = $line =~ /^(.*?).exon./;
      my($linerest) = $line =~ /\t(.*)/;
      my@parents = split("\t",$linerest);
      my$tmpline = "$gene_id";
      push(@pgid, $gene_id);
      print RSLIKE $gene_id;
      print PGID "$gene_id\n";
      foreach (@parents){
        $tmpline .= "\t$_";
        print RSLIKE "\t$_";
      }
      print RSLIKE "\n";
      push(@estgen_input, $tmpline);
    }
  }
  close PGID;
  close RSLIKE;

  my$maskround = "0";

# now start repetitive masking
  while(1){

  $maskround++;
  push(@comments, "starting mask run $maskround in round $round\n");

# retrieve the genomic sequences for the pseudogenes
  foreach my$id(@pgid){
    my$genome_seq_file="$tmpdir/$chr.$id.genome.fa";
    system "/bio/bin/gtf2fa.pl -tx_id $id -allcap -flank 300 $gff_file $master_seq > $genome_seq_file";
  }

# run est_genome reverse: use the parent mRNA vs the Twinscan pseudogene genome sequence
    my@estgen_output = run_alignment(\@estgen_input, $parent, "reverse");
    @estgen_input = @estgen_output;

# with the N.out2 alignment files as input, mask the (masked) genome sequence
    my($false_pos, $flag) = mask_genoseq_allinputs($species);

    if($flag){
# change the master_seq here if it wasn't already done
      $something_was_masked = "1";
      unless($opt_m){
        $master_seq = "$returndir\/$frag.seq.pmasked.masked";
      }
    }
                                                                                                              
# list false positives in the pseudogene file (so we know what happened to the ones
# we do not find in the masked-gtffile)
    if ($false_pos){
      push(@comments, "false positive $false_pos\n");
      open (PSEUDO, ">>$returndir/$frag.$round.pseudo") || graceful_die ("cannot open $pgene_id_list: $!");
      print PSEUDO "\n$false_pos\n";
      close PSEUDO;


# re-create the input files so we're not repetitively trying to mask false positives
      my@false_pos = split("\n", $false_pos);
      foreach(@false_pos){
        my@hold = ();
        my($hit, $gene) = $_ =~ /^(.*?) vs (.*?)/;

# remove gi| etc from NM hit (when running with RefSeq database)
        if($hit =~ /\|/){
          $hit=~ s/.*\|//;
        }
        foreach my$line(@estgen_input){
          if ($line =~ /^$gene/){
             $line =~ s/\t$hit//g;
             unless($line eq "$gene"){
               push(@hold, $line);
             }
          }else{
            push(@hold, $line);
          }
        }
        @estgen_input = @hold;
      }
    }

# quit if opt_m is set
    if($opt_m){
      last;
    }


# The flag contains the size of the masked sequence.
# quit subroutine if nothing was masked or the size is really small 
    if($flag){
      unless ($flag > 20){
        last;
      }
    }else{
      last;
    }
  }
  return "$something_was_masked";
}

# fasta_parse reads in a file with fasta sequences
# creates multiple outputfiles in same directory
# each containint one sequence
# called from run_intron_method

sub fasta_parse{
  my($inputfilename, $ext) = @_;
  my($outdir) = $inputfilename =~ /(.*)\//;
                                                                                                                             
  open (INF, $inputfilename) or graceful_die ("cannot open inputfile $inputfilename: $!");
  my $fasta = new FAlite(*INF);
  while(my $entry = $fasta->nextEntry) {
    my$def = $entry->def;
    my$seq = $entry->seq;
                                                                                                                             
# for EST_GENOME, all sequences must be uppercase
    $seq=~ tr/actg/ACTG/;
                                                                                                                             
    my($name) = $def =~ />(.*?) /;
    unless ($name){
      graceful_die ("cannot find name in $def: format problem for fasta_parse\n");
      next;
    }
                                                                                                                             
    my$outfile = "$outdir\/$chr.$name.$ext";
    open (OUTF, ">$outfile")
      or graceful_die ("cannot open outputfile $outfile in fasta_parse: $!");
    print OUTF "$def\n$seq\n";
    close OUTF;
  }
  close INF;
}


# get_protein_sequence gets the protein nt sequence for an input protein from an input file
# it outputs the seqs to $protein_sequences/ so that can be
# checked first (avoids double work in subsequent runs)
# called from mask_until_oblivion

sub get_protein_sequence{
  my($protein, $protchr, $outfile, $dir) = (@_);
  my$filename = "$dir\/$protchr.gtf";
  my$getseqs_program = "/bio/bin/gtf2fa.pl";
  my$genoseq = "$target_seqdir/$protchr.fa";
  my$found = "0";
  $outfile = "$protein_sequences/$protein.fa";

  if (-e $outfile){
    $found = "1";
  }else{
    push(@comments, "perl $getseqs_program -tx -allcap -tx_id $protein $filename $genoseq > $outfile\n");
    my$out = `perl $getseqs_program -tx -allcap -tx_id $protein $filename $genoseq`;
    if ($out){
      open (OUTF, ">$outfile" ) || graceful_die ("cannot open $outfile: $!");
      print OUTF $out;
      close OUTF;
      $found = "1";
    }else{
      push(@comments, "$protein not found in $species chr$chr"); 
    }
  }
  return $found;

}


# compare blastoutputs compares the informant and prediction hits
# and decides which are nonprocessed pseudogenes, processed pseudogenes and possible paralogs
# called from run_synteny_method

sub compare_blastoutputs{
  my($protein_hits, @informanthits) = (@_);
  my @return = ();
  my @info = ();
  my %possible_nonproces = ();
  my %gene_and_protein = ();	# keeps $species proteins and genes
  my %exons_and_protein = ();	# keeps track of separate exons and protein hits
  my %informant_hits = ();		# keeps track of informant hits per gene
  my %array = ();		# array for possible pseudogenes (introns must be tested for size)
  my %oriprot = ();		# array to remember which (predicted) proteins belong to which fragments

# get and store the protein hits for each exon of each gene
# (there is only one protein hit per gene)
  foreach my$line(@$protein_hits){
    my$found_in_informant = "";
    chomp $line;
    my($id, $prot) = $line =~ /^(.*?)\t(.*?)$/;
    my($gene, $exon) = $id =~ /(.*?)\.exon_(.*)/;

# if the bootstrap method is used (hits start with 'chr'), all parents from the same fragment must be
# treated as one gene (because gene models may split real genes) 
    $oriprot{$id} .= $prot;
    if($prot =~ /^chr/){
      $prot =~ s/\.\d{3}$//;
    }

# keep track of how many exons of a gene hit the same protein
    $gene_and_protein{$gene.$prot}++;

# keep track of the protein hit of each exon
    $exons_and_protein{$id} = $prot;

# see if this exon has a hit in informant
    foreach my$hit(@informanthits){
      if($hit =~ /^$id\,/){
        $found_in_informant = "1";
        $informant_hits{$id} = "1";

# store number of exons per gene that have a hit in informant
        $informant_hits{$gene.$prot}++;
        last;
      }
    }
  }

  my@intronset = ();

# Now, for every exon, see if it is a possible pseudogene
  foreach my$entry(sort keys(%exons_and_protein)){
    my$prot=$exons_and_protein{$entry};
    my($gene, $exon) = $entry =~ /(.*?)\.exon_(.*)/;

# if it's found in informant, it's OK
    if(exists($informant_hits{$entry})){
      push(@info, "$entry (for $prot) found in informant");

# if it's not found in informant, but if there are other exon hits with the same protein 
# AND one of these has a informant hit, it's OK too (?)
    }elsif(($gene_and_protein{$gene.$prot} > 1) && (exists($informant_hits{$gene.$prot}))){
       push (@info, "$entry possible paralog of $prot");

# if there are other exon hits and no informant hits, it's a paralog or a nonprocessed pseudogene,
# or two exons mapping very close together. Put in separate array to find out
    }elsif($gene_and_protein{$gene.$prot} > 1){
       push(@intronset, "$entry.$prot");
    }else{
       push(@return, "$entry parent $oriprot{$entry}");
    }
  }

# check if the exons have small introns between them. This implies that the hits follow each other without gaps
# which is usually the case with multiple exons from one processed pseudogene.
# first split input array in one array per gene

  my%keep_ids = ();
  foreach my$line(@intronset){
    my($gene,$exon, $prot) =$line=~ /(.*)\.exon_([\d]+)\.(.*)/;
# it doesn't matter which of the exons is kept, as long as the original protein is traceable
    $keep_ids{$prot.$gene} = "$exon";
# create new array if one doesn't exist
    unless(defined(@{$array{$gene.'SPLIT'.$prot}})){
      my@tmp = ();
      $array{$gene.'SPLIT'.$prot} = \@tmp;
    }
# add line to array
    my@temp = @{$array{$gene.'SPLIT'.$prot}};
    push(@temp, $line);
    $array{$gene.'SPLIT'.$prot} = \@temp;
  }

# now read every array and test intron- and exonsize with subroutine
  foreach my$entry(sort keys(%array)){
    my$ray = $array{$entry};
    my($gene, $prot) = $entry =~ /(.*)SPLIT(.*)/;
    my$pseudogene = check_introns($ray);
    if($pseudogene){
      my$gene_add=$keep_ids{$prot.$gene};
      my$complete_gene="$gene".'.exon_'."$gene_add";
      push(@return, "$gene multi_exon pseudogene of $oriprot{$complete_gene}");
    }else{
      push(@info, "$gene paralog of $prot");
    }
  }

  return (\@return, \@info);
}

# check_introns aligns intron locations in putative pseudogene and parent
# if introns are not in the same positions, the pseudogene is real
# called from compare_blastoutputs

sub check_introns{
  my($ray) = (@_);
  my$previous_exon = "0";
  my@exons = ();		# array to store exon locations
  my@gene = ();		# array to store gene gff
  my%hs = ();

  foreach my$line(@$ray){
    my$orient;
    my$exonline = "";
    my($gene, $exon) =$line=~ /(.*)\.exon_([\d]+)/;
    my($chr, $twin, $cds, $start, $stop, $score, $sign, $rest);
    if(exists($hs{$gene})){
      if($exon-1 != $previous_exon){
        push (@comments, "$line not following previous\n");
      }

# get exon from gff array
      $orient = $hs{$gene};
      if($orient eq '+'){
        $exonline = $gene[$exon-1];
      }else{
        my$exonpos = $exon * -1;
        $exonline = $gene[$exonpos];
      }
      ($chr, $twin, $cds, $start, $stop, $score, $sign, $rest) = split ("\t", $exonline);
      push(@exons, $start, $stop);
    }else{

# get the correct gff entry
      open(INF, "$gff_file") || graceful_die ("cannot open $gff_file in subroutine check_introns: $!");
      while(<INF>){
        if($_ =~ /$gene\"/){
          unless($_ =~ /codon/){
            ($chr, $twin, $cds, $start, $stop, $score, $sign, $rest) = split ("\t", $_);
            $orient = $sign;
            $hs{$gene} = $sign;
            push(@gene, $_);
          }
        }
      }

# get line for current exon
      if($orient eq '+'){
        $exonline = $gene[$exon-1];
      }else{
        my$exonpos = $exon * -1;
        $exonline = $gene[$exonpos];
      }
      ($chr, $twin, $cds, $start, $stop, $score, $sign, $rest) = split ("\t", $exonline);
      push(@exons, $start, $stop);
    }
    $previous_exon = $exon;
  }

# get the intron and exon sizes and sort. Starts and stops of same exons will
# stay together this way

  sub numerically {$a <=> $b}
  @exons = sort numerically @exons;
  my$previous_stop;
  my@exonsizes = ();
  my@intronsizes = ();
  while(@exons){
    my$start = shift@exons;
    my$stop = shift@exons;
    my$diff = $stop - $start;
    push(@exonsizes, $diff);
    if($previous_stop){
      my$intronsize = $start - $previous_stop;
      push(@intronsizes, $intronsize);
    }
    $previous_stop = $stop;
  }
  my$avg_exon = average (@exonsizes);    
  my$avg_intron = average (@intronsizes);    
  if($avg_intron < $avg_exon){
    return "1";
  }
  return "0";
}


# graceful_exit removes all programs from the tmp directory
# and outputs a reason for exiting in the homedirectory
# called from main and run_est2geno

sub graceful_exit{
  my($reason) = @_;
# remove everyting created in the temporary directory
  unless($opt_t){
    unlink glob("$tmpdir\/$chr.$frag.*");
  }
  unless($opt_c || $opt_m){
# create a N.1.gff link to N.seq.gff if nothing exists
    if(( ! -e "$returndir/$frag.masked.gtf" ) && ( -e "$homedir/$frag.seq.gtf")){
      system "ln -s $frag.seq.gtf $homedir/$frag.1.gtf";
    }else{
      unlink "$homedir/$frag.seq.gtf";
    }
  }


# move all $returndir files to $homedir 
# see if there are any files
  opendir (DIR, "$returndir");
  my@returnfiles=readdir(DIR);
  unless(scalar@returnfiles=="2"){
# move only some returndir files to $homedir if opt_m is set
    if($opt_m){
      system "sort -u $returndir/*.masked.txt >> $homedir/$chr.masked.txt";
      system "cp $returndir/*.masked.gtf $homedir/";
    }else{
      system "rcp $returndir/* $homedir/ " || die "cannot copy!";
    }
  }
  closedir DIR;  

# remove all temporary directories
  unless($opt_t){
    unlink glob ("$returndir/*");
    rmdir "$returndir";
    unlink glob ("$tmpdir/*");
    rmdir "$tmpdir";
  }
  open (INF, ">$homedir\/$chr.pgene.reason") || graceful_die ("cannot open reason file in $homedir:$!");
    print INF "@comments";
    print INF "$reason\n";
  close INF;
  exit(0);
}

# graceful_die replaces die. It creates an error message that tells the user
# where to look for the temporary files. This message is created in an $frag.ERROR file
# in the outputdirectory

sub graceful_die{
  my($reason) = @_;

# create a N.1.gff link to N.seq.gff if nothing exists
print "$round\n";
print "$homedir\n";
print "$frag\n";
  if(($round == "1") && ( -e "$homedir/$frag.seq.gtf")){
#    system "ln -s $frag.seq.gtf $homedir/$frag.1.gtf";
  }else{
#    unlink "$homedir/$frag.seq.gtf";
  }
# move all $returndir files to $homedir using rcp 
#(is nicer than copy using NFS because it waits until the recipient has time)
# see if there are any files
#  opendir (DIR, "$returndir");
# my@returnfiles=readdir(DIR);
#  unless(scalar@returnfiles=="2"){
#    system "rcp $returndir/* $homedir/ " || die "cannot copy!";
#  }
#  closedir DIR;  

# create an error log  
  my $hostname = `hostname`;
  chomp $hostname;
  open (INF, ">$homedir\/$frag.ERROR") || die "cannot open ERROR file in $homedir:$!";
    print INF "@comments";
    print INF "Files can be found in $returndir at $hostname. Please remove files!\n";
    print INF "Died in round $round of $chr\/$frag. Reason: $reason\n";
  close INF;
  exit(0);
}

# get_pseudogenes_from_protein_blastreport fills the coords hash with the 
# genomic protein locations of the parent database
# and uses this hash to identify gene annotations with hits
# elsewhere in the genome (which is indicative of a pseudogene)
# called from run_synteny_method

sub get_pseudogenes_from_protein_blastreport{  
  my($blastfile, $coord, $synteny_file)=(@_);

# the coords hash is used later in the program and has been defined before
  open (COORDS, "$coord") || graceful_die ("cannot open coordinate file $coord: $!");
# some protein have multiple hits
  while(<COORDS>){
    my($id) = $_ =~ /^([\w\.]*?)\t/;
    $coords{$id} .= "NEXT$_";
  }
  close COORDS;
  
  my$qchr = "chr$chr";
  
  open (BLAST, "$blastfile"); 
  my $multiple_report = new BPlite::Multi(\*BLAST);
  
# array for keeping the synteny regions, to be filled with positive hits
  my@all_synteny_hits = ();
  my@hit_ids = ();		# keep track of hit IDs
  my@protein_hits = ();
  my%hs = ();			# keep track of hit gene IDs - add syntenic region
  my%overlap_combi = ();	# keep track of checked overlap combinations

  while(my $blast = $multiple_report->nextReport) {
    my$queryname = $blast->query;
    my($qname) = $queryname =~ /^(.*?) /;
    if($blast->error > 0){
      push(@comments, "No report for $qname\n");
      next;
    }
    (my$frag = $qname) =~ s/\..*//;
    my$samechr_flag = "0";
    my$goodhit_flag = "0";
    my( $q_geno_loc_start, $q_geno_loc_end)= $queryname =~ /start\=(\d*) end\=(\d*)/;
    my$qlength = $blast->queryLength;
    my$secondrun_flag;

    my$subjname;
    while(my $subject = $blast->nextSbjct) {
      my$score;
      my$hsp_length;
      my$percent;
      my$name = $subject->name;

# this must be done first, else the program skips hits
      while( my $hsp = $subject->nextHSP ) {  
        $percent = $hsp->percent;
        $hsp_length = $hsp -> length;
        last;		#only for first HSP
      }

# skip if a RefSeq has a viral description: this is not a pseudogene parent
      if ($name =~ /gag |pol |reverse transcriptase|viral/){
        push(@comments, "WARNING $qname first hit with $name, looks like viral protein, skipping...\n");
        next;
      }
      my$full_name = $name;
      ($subjname) =$name =~ /\>([\w\d\._]*)\s/;

# take only hit that doesn't overlap gene model
# this skips hits with self in bootstrap method
      next if ($subjname =~ /chr$chr\.$frag\./);
# take only hits that are NOT derived from chrUn. Some genes have a double annotation.
      next if ($subjname =~ /^chrUn/);
      next unless ($coords{$subjname});
      my@protein_positions = split("NEXT", $coords{$subjname});
      shift(@protein_positions);
      foreach(@protein_positions){
        chomp $_;
        my$overlap;
        my($name, $subjchr, $subj_geno_start, $subj_geno_end) = split ("\t", $_);

# Sometimes UTR exons overlap parts of the subject that are not positioned on the genome
# therefore, allow a 500k extension on each side of the hit (effectively omitting parents within the same Mb)
        $subj_geno_start -= 500000;
        if($subj_geno_start <1){
          $subj_geno_start="1";
        }
        $subj_geno_end += 500000;

# check overlap of query and subject
        if ($subjchr eq $qchr){
#        if (($subjchr eq $qchr) && ($subjname !~ /^chr[\d\.]+/ )){
# avoid doing extra work
           if(exists($overlap_combi{$subjname.$qname})){
		$overlap="1";
          }else{
            $overlap = check_overlap($frag, $q_geno_loc_start, $q_geno_loc_end, $subj_geno_start, $subj_geno_end);
            $overlap_combi{$subjname.$qname} = $overlap;
          }
          if($overlap){
            $samechr_flag = "1";
          }
        }
      }
# remove really short hits
        if(($percent > 65)&&($hsp_length > 9)){
          $goodhit_flag = "$percent";
        }
# currently only for first hit
      last;
    }
    unless($samechr_flag){
      if($goodhit_flag){
        push(@hit_ids, $qname);
        push(@protein_hits, "$qname\t$subjname");
        (my$gene = $qname) =~ s/\.exon.*//;
# if the gene already had a pseudogene hit, the informant region was already stored
# this is not true for large genes, which may have different informant regions.
# so store separately and remove double hits later
	my@synthits=();
        if(exists($hs{$gene})){
          @synthits=@{$hs{$gene}};
        }
        my$chr_start = $q_geno_loc_start;
        my$chr_end = $q_geno_loc_end;
        unless($opt_m){
          $chr_start += (($frag-1) * 1000000);
          $chr_end +=  (($frag-1) * 1000000);
        }
        my@synteny_hits = is_syntenic($qchr, $chr_start, $chr_end, $synteny_file);
        push(@all_synteny_hits, @synteny_hits);
# get the informant chromosome and MB number per hit
        foreach(@synteny_hits){
          next if ($_ =~ /^$/);
          $_=~ s/chr.*(chr.*\d)/$1/;
          push(@synthits, $_);
        }
        $hs{$gene} = \@synthits;
      }
    }
  }

# remove multiple hits from array
  @all_synteny_hits = reduce_array(@all_synteny_hits); 
  return (\@all_synteny_hits, \@hit_ids, \@protein_hits, \%hs);
}


# create_blastinput selects sequences from an input fasta file
# if their id is present in the input array
# and outputs them in an outputfile
# called from run_synteny_method
  
sub create_blastinput{

  my($protfile, $outf, @hit_ids) = @_;

  open (OUTF, ">$outf") || graceful_die ("cannot open $outf: $!");
  open (PROT, "$protfile") || graceful_die ("cannot open $protfile: $!");
  my $fasta = new FAlite(*PROT);
  while(my $entry = $fasta->nextEntry) {
    my$def=$entry->def;
#      print "$def\n";

    foreach my$id(@hit_ids){
      if($def =~ /^>$id/){
        my$seq=$entry->seq;
        print OUTF "$def\n$seq\n";
        last;
      }
    }
  }
}


# check_overlap takes in a fragment number (in the default setting)
# and two sets of starts and stops. It calculates the genomic position
# from the first set of input coordinates based on the fragment
# then checks if both sets overlap. 
# called from get_pseudogenes_from_protein_blastreport

sub check_overlap{
  my($frag, $twinstart, $twinstop, $RSstart, $RSstop)= @_;
  unless($opt_m){
    my$add = ($frag * 1000000 ) - 1000000;
    $twinstart += $add;
    $twinstop += $add;
  }
  my$overlap = "0";

# see if there is any overlap

# RS             *-----------
# Twin       -----------(-------)

  if(($RSstart >= $twinstart) && ($RSstart <= $twinstop)){
    $overlap = "1";

# RS         ------------*
# Twin    (-------)-----------

  }elsif(($RSstop <= $twinstop) && ($RSstop >= $twinstart)){
    $overlap = "1";

# RS         *--------------*
# Twin          -----------

  }elsif(($RSstart <= $twinstart) && ($RSstop >=$twinstop)){
    $overlap = "1";

# RS         *--------*
# Twin     ---------------

  }elsif(($RSstart >= $twinstart) && ($RSstop <=$twinstop)){
    $overlap = "1";
  }
  return $overlap;

}

# reduce_array removes duplicate enries from an array
# called from get_pseudogenes_from_protein_blastreport

sub reduce_array{
  my@overlap = @_;
  my%hs = ();
  my@single_entries = ();

  foreach(@overlap){
    $hs{$_} = "1";
  }
  foreach (sort keys (%hs)){
    next if ($_=~ /^$/);
    push(@single_entries, $_);
  }
  return @single_entries;
}


# get_good_informant_hits parses blast output for protein vs informant hits
# it selects hits that are from a different part of the genome than the query
# called from run_synteny_method

sub get_good_informant_hits{

  my($blastfile, $syntenic_region) = (@_);
  my@return = ();

  open (BLAST, "$blastfile") || graceful_die ("cannot open $blastfile: $!"); 
  my $multiple_report = new BPlite::Multi(\*BLAST);

  while(my $blast = $multiple_report->nextReport) {
    my$queryname = $blast->query;
    my($qname) = $queryname =~ /^(.*?) /;
    if($blast->error > 0){
      push(@comments, "No report for $qname in good_informant_hits \n");
      next;
    }
    (my$frag = $qname) =~ s/\..*//;
    (my$genename) = $queryname =~ /^(.*?)\.exon/;
    while(my $subject = $blast->nextSbjct) {
      my$score;
      my$hsp_length;
      my$percent;
      my$chance;
      my $name = $subject->name;

# check if hit is with syntenic region that is syntenic with that of the gene
# (mostly important for doing whole chromsome gtfs)
      my$samefrag;
      my ($schr, $syfrag)= $name =~ /\>(.*?)[_\s](\d+)/;
      my$sfrag =(int $syfrag/1000000);
# subjectname contains fragment info. Compare fragment to gene info
      foreach my$entry(@{$$syntenic_region{$genename}}){
	my($syntchr, $syntfrag) = split("\t",$entry);
	if(!$syntchr){
		graceful_die ("no syntchr derivable from $entry, working on $name $qname");
	}
	if(!$schr){
		graceful_die ("no chr found in $name to compare with $entry");
	}
        $syntchr=~s/chr//;
        next unless ($schr eq $syntchr);
        my$mbsyntfrag=(int $syntfrag/1000000);
        if($mbsyntfrag == $sfrag){
           $samefrag="1";
           last;
        }
      }
      while( my $hsp = $subject->nextHSP ) {  
        $percent = $hsp->percent;
        $score = $hsp->score;
        $chance =$hsp->P;
# only take first HSP
        last;
      }
# only take first hit unless it matches to the wrong region
      if($samefrag){
        push(@return, "$qname, $name, $percent percent, score=$score, P=$chance\n");
        last;
      }else{
        push(@comments, "Informant hit $name not in correct syntenic region for $genename\n");
      }
    }	
  }
  return @return;
}


# is_syntenic expects a file with nonoverlapping blocks of synteny
# it outputs an array containing:
# 1 block if the qname region is fully within a synteny block
# if the qname region extends in a non-syntenic region, the previous and/or
# subsequent block is also selected.
# called from get_pseudogenes_from_protein_blastreport

sub is_syntenic{
  my($chr, $qname_geno_start, $qname_geno_end, $synteny_file) = @_;
  my$previous_line="";
  my$hit_line="";
  my$last_line;
  my@informant_hits = ();
  my$found_in_syntenic_region;
  my$add_one_more = "";
# get the corresponding lines from the synteny file
  open (SYN, "$synteny_file" ) ||
    graceful_die ("cannot open synteny file $synteny_file: $!");
  while(<SYN>){
    my$synline=$_;
    next if ($synline =~ /^\#/);
    next if ($synline =~ /chrM/);	# currently skip mitochondrion hits
    my@line = split ("\t", $synline);
    my$humstart = $line[1];
    my$humend = $line[2];

# if we have reached the start but not the end, keep adding
    if($add_one_more){
      push(@informant_hits, $synline);
      if($humend > $qname_geno_end){
        last;
      }
    }
# continue if we're not near the region (keep line)
    elsif($humend < $qname_geno_start){
      $previous_line = $synline;
      next;
# if we pass the start position, add the previous hit
    }elsif($humstart > $qname_geno_start){
      push(@informant_hits, "$previous_line","$synline");
# if we passed the end position, add the next hit
      if($humend >= $qname_geno_end){
        last;
      }else{
        $add_one_more = "1";
        next;
      }
# if we hit the start position, see if the full $species region maps to 1 informant block
    }else{
      if($humstart <= $qname_geno_start && $humend >= $qname_geno_end){
        push(@informant_hits, $synline);
        last;
# if we have the start but not the end, continue until we get the end
      }elsif($humstart > $qname_geno_end){
        push(@informant_hits, $synline);
        last;
      }else{
        push(@informant_hits, $synline);
        $add_one_more = "1";
      }
    }
  }
  close SYN;
  return @informant_hits;
}

# get_informant_coords gets the syntenic chromsome coordinates for putative pseudogenes
# called from run_synteny_method

sub get_informant_coords{
  my($informant_hits) = (@_);

# some informant fragments may be adjacent. Create one big fragment out of those
  my$previous_end = "0";
  my$previous_start = "0";
  my$previous_chrnr = "0";
  my$newstart = "0";
  my$newend = "0";
  my$chrnr = "";
  my@informant_coords = ();

  foreach(@$informant_hits){
    chomp $_;
    my@frags_between;	# array for keeping track of complete fragments in synteny region
    my@line = split ("\t", $_);
    my$musstart = $line[4];
    my$musend = $line[5];
    ($chrnr = $line[3]) =~ s/chr//;
    if(($musstart == $previous_end) && ($chrnr eq $previous_chrnr)){
      $newend = $musend;
      $newstart =$previous_start;
      $previous_end = $musend;
    }else{

      if($newend){
        push (@informant_coords,  "$previous_chrnr\t$newstart\t$newend");
      }
      $previous_chrnr = $chrnr;
      $previous_end = $musend;
      $previous_start = $musstart;
      $newend = $musend;
      $newstart = $musstart;
    }
  }
# print last line
  push (@informant_coords,  "$previous_chrnr\t$newstart\t$newend");
  return @informant_coords;
}

# get_intronsizes calculates intronsizes from an array with exonlocations
# called from mask_genoseq_allinputs

sub get_intronsizes{
  my($exons) = (@_);
  my$previous_end;
  my(@introns) = ();

  foreach my$start(sort numerically keys(%$exons)){
    my$end = $$exons{$start};
    if($previous_end){
      my$intronsize = $start - $previous_end;
# test for bad alignment
      if($intronsize < "0"){
        @introns = ();
        push(@introns, "negative intronsize");
        return @introns;
      }
      push(@introns, $intronsize);
    }
    $previous_end = $end;
  }
  return @introns;
}

# average averages the numbers in an array
# called from mask_genoseq_allinputs and check_introns

sub average{
  my@numbers = @_;
  my$nr = scalar @numbers;
  unless($nr){
    $nr = "1";
  }
  my$total = "0";
  foreach(@numbers){
    $total += $_;
  }
  my$avg = sprintf('%.2f', $total/$nr);
  return $avg;
}


# run_intron_method is the intron alignment method of pseudogene finding
# called from main

sub run_intron_method{

  while (1){

  push(@comments, "starting round $round using intron alignment method\n");

# create intron coordinates file
    my$intronfile = "$returndir\/$frag.$round.introns";

    unless (-e "$intronfile"){
      my$introns = intron_coords("$gff_file");
      open (INTRONS, ">$intronfile") || graceful_die ("cannot create $intronfile: $!");
      print INTRONS "$introns";
      close INTRONS;
    }

# get gene sequences and store in $tmpdir
    my$geneseqs="$tmpdir\/$chr.$frag.genes.fa";
    if (-e "$tmpdir/$frag.seq.masked"){
      system "/bio/bin/gtf2fa.pl -allcap -tx -cds $gff_file $tmpdir/$frag.seq.masked > $geneseqs";
    }else{
      graceful_die ("cannot find $tmpdir\/$frag.seq.masked");
    }

# get the genomic locations from this file when running whole chr files
    if($opt_m){
      my@geno_from_file=`grep '>' $geneseqs`;
      foreach(@geno_from_file){
        my($id, $start)=$_=~/>(.*?) start\:(\d+) end/;
	my$clonefrag= int ($start/1000000)+1;
	$clone_frag{$id}=$clonefrag;
      }
    }

# unsplit the file (gene sequences are put in same dir as $geneseqs file)
    fasta_parse($geneseqs, "fa");

# run Blast
    my$blast_flags = "-E=0.0001 -M=1 -N=-1 Q=2 R=2 -cpus=1 -warnings -notes";
    my$blastout = "$tmpdir/$chr.$frag.blastout";

    push(@comments, "running blast using $parent_blastfile, output in $blastout\n");
    system "$BLASTN $parent_blastfile $geneseqs $blast_flags > $blastout";

# remove inputfile
    unless($opt_t){
      unlink $geneseqs;
      push(@comments, "$geneseqs removed\n");
    }

# parse blast output
    my$rs_hitfile = "$returndir\/$frag.$round.RS_hits";
    my($rs_hits, $ng_hits)=get_parent_hits_from_blast ($blastout, $round, $species);
    my($ng_maskflag, $rs_maskflag);

# remove the blastoutput file: it can be very big
    unless($opt_t){
      unlink "$blastout";
    }

# quit if there are no hits in the blastfile
    unless($$rs_hits[0]){
      unless($opt_t){
        unlink glob("$tmpdir\/$chr.$frag.*");
      }
      return("no parent blast hits in round $round");    
    }
# define outputfile
    my$pseudogenes_file="$returndir\/$frag.$round.pseudo";

# compare introns if normal parents hits are found
    if($$rs_hits[0]){
      my@newintrons;
      open (RSOUT, ">$rs_hitfile");
      foreach(@$rs_hits){
        print RSOUT "$_";
      }
      close RSOUT;

# run est2genome for each gene-parent combination 
# change alignment program (because it may be sim4)
# output array is not important here
      my$hold=$alignment_program;
      $alignment_program = "est2geno";
      my@obsolete = run_alignment($rs_hits, "refseq", "forward");
      $alignment_program = $hold;

# compare the introns in the est_genome outputfiles to the intron positions
# in N.$round.introns
      @newintrons = intron_diffs_in_est2geno($intronfile);

# mask if there are pseudogene hits
      if($newintrons[0]){
        open (PSEUDO, ">>$pseudogenes_file");
        print PSEUDO @newintrons;
        close PSEUDO;
        $rs_maskflag=mask_until_oblivion(\@newintrons, "refseq");
      }else{

# quit if there are no pseudogene hits after intron finding
        unless($opt_t){
          unlink glob("$tmpdir\/$chr.$frag.*");
        }
        return("no additional introns: no pseudogenes in round $round");
      }
    }

# quit if nothing was masked
    unless($rs_maskflag){
      return("No additional masking in round $round");
    }

  if($opt_n){
    graceful_exit("Option n set: exiting before Twinscan run");
  }elsif($opt_m){
    $round++;
    return("Option m set");
  }

# with the newly masked file, run Twinscan
  run_twinscan();

# add one to counter, and rerun
    $round++;
    if ($round > $round_cutoff){
      graceful_exit("ERROR: round number too high!");
    }
  }
}

# get_parent_hits_from_blast parses blast output and
# retrieves 'good' hits. Outputs an array containing NM hits
# called from run_intron_method

sub get_parent_hits_from_blast { 

  my($infile, $round, $species) = @_;
  
  my @introns = ();
  my @ng_hits = ();
  open (BLAST, "$infile"); 
  my $multiple_report = new BPlite::Multi(\*BLAST);
  my $nr_of_genes = 0;
  my $double_orientation = 0;

  while(my $blast = $multiple_report->nextReport) {
    my$queryname = $blast->query;
    my($qname) = $queryname =~ /^(.*?) /;
    if($blast->error > 0){
      push(@comments, "No report for $qname\n");
      next;
    }
    my$qlength = $blast->queryLength;
    $nr_of_genes++;
    my $nr_of_hits = "0";
    my$maxscore;
    my $i = "0";
    my $ps_flag = "0";
    my $matches = "";
    my $ng_matches = "";
    my @covered = ();
    my @oldcovered = ();
    my $cutoff = "2";

# for whole chromosome input from external sources, fragment cannot be derived from gene name
# get fragment from hash to compare with chr id in bootstrap method. 
    my$clonefrag=$clone_frag{$qname};

    while(my $subject = $blast->nextSbjct) {
      my$score;
      my$percent;
      my $nr_of_pseudogenes = 0;
      $i++;
      $nr_of_hits++;
      my@onehit_covered = ();
      my $name = $subject->name;
      my$full_name = $name;
# catch the NM_xxx ID
      my($subjname) = $name =~ /(N\w_[\d\.]*)/;
# allow for non-RefSeq IDs
      unless($subjname){
        ($subjname) = $name =~ />(.*?) /;
      }

      my$percent_flag = "0";
      my$minusstrand = "0";

      my@covered_region = ();
      while( my $hsp = $subject->nextHSP ) {
# take only score for first HSP
        unless($score){
          $score = $hsp->score;
        }
        my$qstart=$hsp->queryBegin;
        my$qend=$hsp->queryEnd;
        push(@covered_region, $qstart, $qend);

        my$percent = $hsp->percent;
        if($percent > 75){
          $percent_flag = "1";
        }
# skip all hits that are Minus / Plus (these are genes on the other strand
        if($qend < $qstart){
          $minusstrand = 1;
        }
      }

      next unless($percent_flag);
      next if ($minusstrand);

      if($opt_m){
        next if ($subjname =~ /chr$chr\.$clonefrag\./);
      }

# set maxscore with first good hit found
      unless($maxscore){
        $maxscore = $score;
      }
      push(@onehit_covered, @covered_region);


# only work with NM hits when using RefSeq
      if ($subjname =~ /N[GRC]_/){
        next;
      }else{

# create a cutoff (low scores should not be reported)
        if ($score/$maxscore > "0.75"){
          $matches .= "\t$subjname";
        }
# First, create a nonoverlapping range set for all HSPs 
# found with this subject
        my@reduce = reduce_ranges(@onehit_covered);
# Add the range set to the existing ranges...
        @oldcovered = @covered;
        push(@covered, @reduce);
      }
# ... and check again
      my@reduce = reduce_ranges(@covered);
      my@reducer = reduce_ranges(@reduce);
      @covered = @reducer;
      my$size = scalar(@reducer);

# check if the new subject has a different range
      if((scalar(@reducer) > scalar(@oldcovered))&&($i>1)){
         unless($matches =~ /$subjname/){
           $matches.= "\t$subjname";
         }
        my$refscore = $score/$maxscore;
        my$rounded = sprintf("%.2f", $refscore);
        $cutoff += "2";
      }
      if((scalar(@reducer) < scalar(@oldcovered))&&($i>1)){
        next;
      }

# check if the new subject extended any range by more than 100 nt
      while(@oldcovered){
        my$oldstart = shift@oldcovered;
        my$oldend = shift@oldcovered;
        my$newstart = shift@reducer;
        my$newend = shift@reducer;
        next if(($newstart == $oldstart)&&($newend == $oldend));
        if (($oldstart-$newstart) > 100){
          unless($matches =~ /$subjname/){
            $matches.= "\t$subjname";

          }
        }
        if(($newend - $oldend) > 100){
         unless($matches =~ /$subjname/){
           $matches.= "\t$subjname";
         }
        }
      }

      while( my $hsp = $subject->nextHSP ) {  } # if this is not in, the program fails
    }
    if($matches){
      my$line="$qname"."$matches"."\n";
      push (@introns, $line);
    }elsif($ng_matches){
# output NG matches only if there is no good parent gene
      my$line="$qname"."$ng_matches"."\n";
      push (@ng_hits, $line);
    }else{
    }
  }
  close BLAST;
  return (\@introns, \@ng_hits);
}


# intron_diffs_in_est2geno takes in an inputfile containing gene predictions 
# and their intron positions, and a directory containing est_genome outputfiles
# it outputs a list of gene prediction IDs that show extra introns in the
# alignments
# called from run_intron_method

sub intron_diffs_in_est2geno{
  my($coordsfile) = @_;
  my$est;
  my$genome;
  my$nr_of_pred_exons;
  my$old_est = "";
  my$parents = "";	
  my$splice_flag;
  my@allparents = ();
  my@files = ();
  my@return = ();
  my%hitexons = ();
  my%list_of_genes= (); 

# create table with intronsizes (from Twinscan gtf file)
  open (COORDS, "$coordsfile") ||
    graceful_die ("cannot open file $coordsfile in intron_diffs_in_est2geno: $!"); 
  my@sizelist=<COORDS>;
  close COORDS;

# must sort files before reading them
  opendir (DIR, "$tmpdir");
  while (defined(my$file = readdir(DIR))){
    if(($file eq '.') || ($file eq '..')){next}
    unless($file =~ /^$chr.*\.out$/){next}
    push(@files, $file); 
  }
  closedir DIR;  
  
  my@file = sort(@files);  
  foreach my$file(@file){ 
    open (INF, "$tmpdir\/$file")|| graceful_die ("cannot open $file: $!"); 
    if(my $report = new EGlite(\*INF)){
      unless($report->hasIntron) {next}
      $est= $report->est;
      $genome= $report->genome;
      $genome =~ s/\.[0-9]$//;	# remove version number

# see if the report has no -Introns
# and store the ?Introns
      my$intron_count="0";
      my$check_boundaries;
      my%boundary=();
  
      while(my $element = $report->nextElement) {
        if($element->type eq "INTRON"){
          $intron_count++;
          my$tmpdirect = $element->direction;
  
# if all boundaries are bad, do not continue
          if($tmpdirect eq "+"){
            $check_boundaries = "continue";
  
# if there is a -Intron, skip the alignment altogether
          }elsif($tmpdirect eq "-"){
            push(@comments, "minus intron in $est vs $genome, skipping\n");
            last;
          }else{
            $boundary{$intron_count} = "1";
          }
        }
      }
      next unless($check_boundaries);

      my$estline = $est;
      my$predicted_introns;
      my%pred_int = ();
      my@pred_int = ();
  
# see if we're already working on this est
      if ($est eq $old_est){
# skip to last part
      }else{
  
# print out results of former est
        if($splice_flag){
          my@exons = keys(%hitexons);
            foreach my$ex(@exons){
              push(@return, "$old_est\.exon_$ex\t$hitexons{$ex}\n");
              push(@allparents, "$hitexons{$ex} ");
          }
  
# the last exon does not have a boundary, so:
          unless(@exons){
            push(@return, "$old_est\.exon_$nr_of_pred_exons\t$parents\n");
              push(@allparents, "$parents ");
          }
        }
  
# re-initialize 
        $splice_flag = "0";
        %hitexons = ();
        $parents = "";
      }

# find corresponding intron sizes for est
      foreach(@sizelist){
        if ($_ =~ /^$est/){
          $predicted_introns = $_;
          chomp$predicted_introns;
          @pred_int = split("\t", $predicted_introns);
  
# remove the gene name
          shift(@pred_int);
          $nr_of_pred_exons = scalar(@pred_int)+1;
          my$i = "0";
          foreach(@pred_int){
            $i++;
            $pred_int{$_} = $i;
          }
        }
        if($predicted_introns){last}
      }
      unless($predicted_introns){
        push (@comments, "no intron list found for $est, skipping...\n");
        next;
      }
      $old_est = $est;
      
      my $est_begin = $report->estBegin;
      my $est_end = $report->estEnd;
      my $exon_nr = "0";
      my $nr_of_exons = "0";
  
# counter to make sure we skip the last exon, because the last
# estEnd is not the start of an intron
      my$nr_of_elements = $report->countElements;
  
      $report->resetElements;
      while ($report->hasMoreElements){
      my $element = $report->nextElement; 
        $nr_of_elements--;
        if (($element->type eq "EXON")&&($nr_of_elements)){
          $exon_nr++;
  
# skip exons with a bad intron/exon boundary
          if(exists($boundary{$exon_nr})){
            next;
          }
# skip exons with a bad or negative score
          my$exon_score = $element->score;
          if($exon_score < 0){
            next;
          }
          my$exon_length = $element->length;
          my$div = $exon_score/$exon_length;
          if($div < "0.40"){
            next;
          }
          my$intron_loc = $element->estEnd;
  
# see if the intron/exon boundary exists
          if(exists($pred_int{$intron_loc})){
          }else{
            $splice_flag = "1";
            unless(@pred_int){
              $hitexons{single_exon_gene} = "$genome";
            }
            foreach(@pred_int){
              if($_ > $intron_loc){
                my$exon_nr = $pred_int{$_};
                if(exists($hitexons{$exon_nr})){
                  unless($hitexons{$exon_nr} =~ /$genome/){
                    $hitexons{$exon_nr} .= " $genome";
                  }
                }else{
                  $hitexons{$exon_nr} = "$genome";
                } 
                last;
              }
            }
# keep track of parent hits
            unless($parents =~ /$genome/){
              $parents .= "$genome ";
            }
          }
        }
      }
    }
    close INF; 
  } 
  
# now print out results of last est
  if($splice_flag){
    my@exons = keys(%hitexons);
    foreach my$ex(@exons){
      push(@return, "$old_est\.exon_$ex\t$hitexons{$ex}\n");
      push(@allparents, "$hitexons{$ex} ");
    }
    unless(@exons){
      push (@return, "$old_est\.exon_$nr_of_pred_exons\t$parents\n");
      push(@allparents, "$parents ");
    }
  }
  return(@return);
}

# intron_coords gets the intron coordinates in a gene
# based on a gtf file
# called from run_intron_method

sub intron_coords{
  my ($filename)=@_;
  my $LINE = 120;
  my $gtf = GTF::new({gtf_filename => $filename});
  my $return = "";
		    
  my $genes = $gtf->transcripts;
  foreach my $gene (@$genes){
    my$tmpid=$gene->id;
#start is -1 if there is no start
    next if ($gene -> start < "1");
    my $exons = $gene->cds;
    my $minus_strand = 0;
    $return.= $tmpid;
    if($gene->strand eq "-"){
      next unless (exists $$exons[$#$exons]);
      my $pos = $$exons[$#$exons]->length;
      for(my $i = $#$exons - 1; $i >= 0;$i--){
        $return .= "\t". $pos; 
        $pos += $$exons[$i]->length;
      }
    }
    else{
        next unless (exists $$exons[0]);
	my $pos = $$exons[0]->length;
	for(my $i = 1; $i <= $#$exons;$i++){
          $return .= "\t". $pos; 
          $pos += $$exons[$i]->length;
        }
     }	
     $return.= "\n";
  }
  return $return;
}

# reduce_ranges takes in an array containing number pairs
# the number pairs are starts and ends of hits. The routine looks for
# overlapping hits and merges these
# called from get_parent_hits_from_blast
sub reduce_ranges{

  my@ranges = @_;
  my@reduced = ();
  my$length = scalar(@ranges);
  my$count = "0";
  my%pairs = ();
  my$warning = " ";

  while($count<$length){

    my$newstart = $ranges[$count];
    my$newend = $ranges[($count+1)];
    $count+="2";
    my$new_pair = "0";
    unless(%pairs){
      $pairs{1}{start}=$newstart;
      $pairs{1}{end}=$newend;
      next;
    }
    my@keys = keys %pairs;
    my$range_nr = scalar@keys;
    foreach my$key(@keys){
      $new_pair = "0";
      my$start = $pairs{$key}{start};
      my$end = $pairs{$key}{end};
      if(($newstart<=$end)&&($newstart>=$start)){	# start lies in range
        if($newend > $end){			# end extends range
          $count = "0";
# check if the extension is significant
          $pairs{$key}{end}= $newend;
        }else{last}
      }elsif(($newend<=$end)&&($newend>=$start)){	# end lies in range
        if($newstart < $start){			# start extends range
          $count = "0";
          $pairs{$key}{start}= $newstart;
        }else{last}
      }elsif(($newend>$end)&&($newstart<$start)){	# both extend range
          $count = "0";
          $pairs{$key}{start}= $newstart;
          $pairs{$key}{end}= $newend;
          last;
      }elsif(($newend<=$end)&&($newstart>=$start)){	# both within range
	last;						# do nothing
      }else{					# range does not overlap		
        $new_pair = 1;    
      }
    }

    if($new_pair){
      $pairs{($range_nr+1)}{start} = $newstart;
      $pairs{($range_nr+1)}{end} = $newend;
    }

  }
  my@keys = keys %pairs;
  foreach my$key(@keys){
    push(@reduced, $pairs{$key}{start}, $pairs{$key}{end});
  }
  return @reduced;
}

# run_alignment runs est_genome for all sequences in a list
# alignment program is est2geno or sim4 (selectable)
# outputfiles are put in temporary directory
# returns a reduced list with only the sequences that actually aligned.
# called from mask_until_oblivion

sub run_alignment{
  my($inf, $parent, $orient) = (@_);
  my@return = ();
  my%hs = ();  
                                                                                                                             
  foreach my$line(@$inf){
# split the line in gene and parent IDs
    my(@seqnames) =  split (" ", $line);
                                                                                                                             
# the first entry is the gene ID
    my$twin = shift@seqnames;

# the rest of the entries are parents: put in hash
     while(@seqnames){
      my$rs_id = shift@seqnames;
      $hs{$twin}{$rs_id} = "1";
#      print "seqname is $rs_id\n";
    }
  }

# now for every prediction, take all the parents, check existence and compare sizes
  foreach my$twin(sort keys %hs){

# set the maximum rs size 
    my$max_rs_size = "0";
    foreach my$rs_id(sort keys %{$hs{$twin}}){
      my$rs_seq = "$parent_geno_seqs/$rs_id\.fa";

# only care about size if it's a mapback
      if($orient eq "reverse"){
        $rs_seq = "$parent_clone_seqs/$rs_id\.fa";
        if(($parent eq "protein") && ($rs_id !~ /^chr/)){
          $rs_seq = "$protein_sequences/$rs_id.fa";
        }
        if ( -e "$rs_seq" ){
# list size of parent and keep if it's larger than the maximum parentsize found for this gene model
          open(EST, "$rs_seq");
          my $fasta = new FAlite(*EST);
          while(my $entry = $fasta->nextEntry) {
            my$seq = $entry->seq;
            my$estsize = length$seq;
            if($estsize > $max_rs_size){
              $max_rs_size = $estsize;
            }
          }

# rs sequence not found, skip
        }else{
          push (@comments, "could not find $rs_seq for alignment, skipping...\n");
          delete $hs{$twin}{$rs_id};
        }
      }
    }

# continue only if a max_rs_size was set (only happens in "reverse")
    if($max_rs_size){

      my$twinseq = "$tmpdir\/$chr.$twin\.fa";
      if($orient eq "reverse"){
        $twinseq = "$tmpdir\/$chr.$twin\.genome.fa";
      }

# get size of prediction
      my$genosize = "0";
      open(GEN, "$twinseq") || graceful_die ("cannot open $twinseq in mask subroutine: $!");
      my $fasta2 = new FAlite(*GEN);
      while(my $entry = $fasta2->nextEntry) {
        my$seq = $entry->seq;
        $genosize = length$seq;
      }
# if the prediction sequence is too small, enlarge
      if(3*$max_rs_size > $genosize){
                                                                                                                             
# enlarge Twinscan genome file
        my$enlarge_id = "$tmpdir\/$chr.$frag.enlarge.txt";
        open (TMPF, ">$enlarge_id");
        print TMPF "$twin\n";
        close TMPF;
        my$genome_seq_file="$tmpdir/$chr.$twin.genome.fa";
        my$mastergff = "$returndir/$frag.seq.gtf";
        system "/bio/bin/gtf2fa.pl -tx_id $twin -allcap -flank $max_rs_size $mastergff $master_seq > $genome_seq_file\n";
        push(@comments, "command is /bio/bin/gtf2fa.pl -tx_id $twin -allcap -flank $max_rs_size $mastergff $master_seq > $genome_seq_file\n");

        unlink "$enlarge_id";
      }
    }
  }

# now the file size should be OK and all rs files needed should be there.
# run alignment for all sequences remaining in the hash
# First, establish program. (Set as 'our' variable in main program)

# space parameter is set to allow for bigger memory use - this should use max 1 Gig
  my$est2geno_params = "-splice_penalty 10 -intron_penalty 20 -mismatch 3 -space 1000 -minscore 10 -align";

  foreach my$twin(sort keys %hs){
    foreach my$rs_id(sort keys %{$hs{$twin}}){

      my$twinseq = "$tmpdir\/$chr.$twin\.fa";
      my$checkfile = "$tmpdir\/$chr.$twin"."_"."$rs_id.out";
      my$rs_seq = "$parent_geno_seqs/$rs_id\.fa";
      if($orient eq "reverse"){
        $twinseq = "$tmpdir\/$chr.$twin\.genome.fa";
        $checkfile = "$tmpdir\/$chr.$twin"."_SPLIT_"."$rs_id.out2";
        $rs_seq = "$parent_clone_seqs/$rs_id\.fa";
        if(($parent eq "protein") && ($rs_id !~ /^chr/)){
          $rs_seq = "$protein_sequences/$rs_id.fa";
        }
                                                                                                                             
        if($alignment_program eq "sim4"){
          system "$SIM4 $rs_seq $twinseq A=4 R=0 1>$checkfile 2>/dev/null";
        }elsif($alignment_program eq "est2geno"){
          system "$EST2GENO -est $rs_seq -genome $twinseq $est2geno_params > $checkfile";
        }
# check if an alignment came out; remove rs from hash if not
        my$got_align;
        open(TMP, "$checkfile");
        while(<TMP>){
          if ($_ =~ /\|\|\|/){
            $got_align = "1";
          }
        }
        unless($got_align){
          delete $hs{$twin}{$rs_id};
          push(@comments, "WARNING $checkfile came out empty in round $round\n");
          unlink ("$checkfile")|| graceful_die ("cannot remove $checkfile: $!");
        }

      }else{
# orient is "forward", runs only est2geno MUST BE CHANGED
        system "$EST2GENO -est $twinseq -genome $rs_seq $est2geno_params > $checkfile";
      }
    }
  }
# of the remaining hash entries, create return array
  foreach my$twin(sort keys %hs){
    my$line = "$twin";
    foreach my$rs(sort keys %{$hs{$twin}}){
      $line .= "\t$rs";
    }
    push(@return, $line);
  }
  return @return;
}

# mask_genoseq_allinputs reads est2geno out2 files or sim4 out2 files based on input
# and masks everyting found in the file. 
# called from mask_until_oblivion

sub mask_genoseq_allinputs{
  my($species) = @_;

  my$parser = "EGlite";
  if($alignment_program eq "sim4"){
    $parser = "Sim4Parser";
  }

  push(@comments, "using $parser for masking\n");

# this outputfile will be created if it doesn't exist. If it does exist, $master_seq is the same as this outputfile
  my$maskedfile = "$returndir\/$frag\.seq\.pmasked\.masked";
  my$maskgtf = "$returndir\/$frag.masked.gtf";
  my$masktext = "$returndir/$frag.$round.masked.txt";

  my$false_pos = "";	#will be returned to main program
  my$mask_flag = "0";	#changes to 1 as soon as something is masked
  my@mask_text = ();
  my@mask_gtf = ();

# read in all *.out2 files
  opendir (DIR, "$tmpdir");
  my@allfiles = grep {/out2$/} readdir(DIR);

# reverse sort so the exon alignment files get read last
  @allfiles = reverse(@allfiles);
  foreach my$file(@allfiles){

# some files may be deleted througout the subroutine
    next unless ( -e "$tmpdir/$file");
    open (INF, "$tmpdir\/$file")|| graceful_die ("cannot open $file: $!");
    push(@comments, "working on $tmpdir\/$file\n");
    
    if(my $report = new $parser(\*INF)){
      my$genome = $report->genome;
      $genome =~ s/\s.*//;
      my$cDNA = $report->est;
      my($est) = $cDNA =~ /^(.*?)\s/;
# allow for proteins
      unless($est){
        $est = $cDNA;
      }
      $est =~ s/\|$//;
      if($est =~/\|/){
        $est =~ s/.*\|//;
      }
      (my$checkfile = $file) =~ s/_SPLIT_(.*)out2$/\.genome\.fa/;
# find start position of genome sequence on chromosome
      unless (open (TWINGEN, "$tmpdir/$checkfile")){
        push (@comments, "can't open $checkfile genome file, skipping...\n");
        next;
      }
      my$line = <TWINGEN>;
      my($loc1, $loc2, $strand) = $line =~ 
        /\>.*?\sstart\:(\d*)\send\:(\d*) strand (.)/;
      graceful_die ("something wrong with first line of file $genome") unless($loc2);
      close TWINGEN;
      my%maskexons = ();
      my@intron_sizes=();
      my@exon_sizes=();
      my@intron_coords = ();

# select the regions that must be masked in genome seq
      while(my $element = $report->nextElement) {

        my($geno_start, $geno_end);
# get intronsize (works only in est2geno)
        if($element->type eq "INTRON"){
          $geno_start = $element->genomeBegin;
          $geno_end = $element->genomeEnd;
          my$intronsize = $geno_end - $geno_start;
          if ($strand eq "-"){
            $intronsize = $geno_start - $geno_end;
          }
          push(@intron_sizes, $intronsize);
          next;
        }

        $geno_start = $element->genomeBegin;
        $geno_end = $element->genomeEnd;
# keep the original exon locations in case we find an interrupted alignment
        push(@intron_coords, $geno_start, $geno_end);

# now get the genomic locations
        my$maskbegin = $loc1+$geno_start;
        my$maskend = $loc1+$geno_end;

# if the Twinscan gene is on the minusstrand, the sequence numbers are 
# reversed and the positions should be subtracted
        if ($strand eq "-"){
          $maskend = $loc2 - $geno_start +2;
          $maskbegin = $loc2 - $geno_end +2;
        }

# store start and end for masking out later
        $maskexons{$maskbegin} = $maskend;
        my$length = $maskend - $maskbegin;
        push(@exon_sizes, $length);
      }
# if @exon_sizes only contains very small fragments, do not mask
      if((scalar@exon_sizes == "1") && ($exon_sizes[0] < 50)){
        $false_pos .= "$est vs $genome: Single hit, size smaller than 50, skipping...\n";
        unlink ("$tmpdir\/$file");
        next;
      }

# intron sizes are not explicity listed in the Sim4 report. Extract them from 
# exonsizes with a subroutine
      if($alignment_program eq "sim4"){
        @intron_sizes = get_intronsizes(\%maskexons);
        if(exists($intron_sizes[0]) && ($intron_sizes[0] =~ /negative/)){
# bad alignment, do not mask
          $false_pos .= "$est vs $genome has bad alignment, skipping...\n";
          unlink ("$tmpdir\/$file");
          next;
        }
      }

# If the alignment has more than one exon, see if it looks like a real gene
# (small exons, large introns). If so, do not mask.
      if(scalar@exon_sizes > 1 ){
        my$total_realgene_flag = "";
        my$total_exon_avg = average (@exon_sizes);
        my$total_intron_avg = average (@intron_sizes);

# chicken (and worm) has smaller introns
        if (($total_exon_avg < $total_intron_avg) && ($species eq "chicken"|| $species eq "worm")){
          push(@comments, "$est vs $genome total exon avg: $total_exon_avg intron avg: $total_intron_avg\n");
# Before marking this a false positive, first see if the 'intron' is not all repeat
# (some pgenes are interrupted by repeats);
          my$nflag = check_for_N(\@intron_coords, "$tmpdir/$checkfile");
          unless($nflag){
            push(@comments, "setting total realgene flag...\n");
            $total_realgene_flag = "1";
#            my($predgene, $prot) = $file =~ /^($chr.$frag.\d{3}).*?_(.*)\.out2/;
            my($predgene, $prot) = $file =~ /^(.*?)_SPLIT_(.*)\.out2/;
            push(@comments, "removing all $predgene $prot (from $file) combination out2 files\n");
            unlink glob ("$tmpdir\/$predgene*$prot.out2");
            next;
          }
        }elsif ((2*$total_exon_avg) < $total_intron_avg){
          push(@comments, "$est vs $genome total exon avg: $total_exon_avg intron avg: $total_intron_avg\n");
# Before marking this a false positive, first see if the 'intron' is not all repeat
# (some pgenes are interrupted by repeats);
          my$nflag = check_for_N(\@intron_coords, "$tmpdir/$checkfile");
          unless($nflag){
            push(@comments, "setting total realgene flag...\n");
            $total_realgene_flag = "1";
            my($predgene, $prot) = $file =~ /^($chr\..*?)_SPLIT_(.*)\.out2/;
            push(@comments, "removing all $predgene $prot ($chr from $file) combination out2 files\n");
            unlink glob ("$tmpdir\/$predgene*$prot.out2");
            next;
          }
        }
      }
# if the routine gets through to here, the est is ok for masking out
# open sequence file for masking
# make sure all masked regions end up in the same file

      (my$hitgene) = $file =~ /^(.*?)_SPLIT_/;

# do not mask with whole chromsome option
        my$name;
        my$seq;
      unless($opt_m){
        open (GEN, "$master_seq") || graceful_die ("cannot open file $master_seq: $!");
        my $fasta = new FAlite(*GEN);
        while(my $entry = $fasta->nextEntry) {
          $name=$entry->def;
          $seq=$entry->seq;
        }
        close GEN;
      }

      my@starts = sort keys(%maskexons);
      foreach my$begin(@starts){
        my$end = $maskexons{$begin};
        my$masklength = $end - $begin;
        $mask_flag = "$masklength";
        my$replace = ("N")x($masklength+1);

# create GTF maskfile and textfile
        push (@mask_gtf, 
          "$frag\tPMASK\tCDS\t$begin\t$end\t\.\t+\t\.\tgene_id \"$frag.$est\"\; transcript_id \"$frag.$est\.1\"\;\n");
        push (@mask_text, "$hitgene\n");
        unless($opt_m){
          substr($seq, $begin-2, $masklength+1) = $replace;
        }
      }
      unless($opt_m){
        open (MASK, ">$maskedfile") || graceful_die ("cannot open seq outputfile: $!");
        print MASK "$name\n$seq\n";
        close MASK;
# change masterseq after first time masking
        $master_seq = "$returndir\/$frag.seq.pmasked.masked";
      }

# if there is no readable file in the *out2: issue warning
    }else{
      push(@comments, "$tmpdir\/$file is empty or has wrong format\n"); 
    }
    close INF;
# remove file if it still exists
    if ( -e "$tmpdir\/$file" ){
      unlink ("$tmpdir\/$file");
    }
  }
  closedir DIR;
# append masked regions to outputfile
  if($mask_gtf[0]){
    open (INF, ">>$maskgtf");
    print INF @mask_gtf;
    close INF;
# append masked regions to outputfile
    open (INF2, ">>$masktext");
    print INF2 @mask_text;
    close INF2;
  }
  return ($false_pos, $mask_flag);
}

# check_for_N is called when a large interruption is found in a
# pseudogene-to-genome alignment. It takes in the genome sequence and an array containing the
# starts and ends of all the exons in the alignment. These are used as ends and starts of
# intron locations. The intron is extracted from the file and the fraction of Ns is determined
# if the fraction is > 75%, the subroutine returns 'true'
# count ALL introns!
# called from mask_genoseq_allinputs
                                                                                                                             
sub check_for_N{
  my$return = "0";
  my($intron_coords, $checkfile) = (@_);
  my$intronseq="";
  my$intronsize="1";

# remove first and last of array (because it starts with an end and ends with a start)
  shift@$intron_coords;
  pop@$intron_coords;
  open (SEQFILE, "$checkfile") || die "cannot open $checkfile: $!";
  my $fasta = new FAlite(*SEQFILE);
  while(my $entry = $fasta->nextEntry) {
    my$seq = $entry->seq;
    while(@$intron_coords){
      my$start = shift@$intron_coords;
      my$end = shift@$intron_coords;
      my$size = $end - $start;
      $intronsize+=$size;
      my$seq = substr($seq, $start, $size-1);
# add sequence to total intron sequence
      $intronseq.=$seq;
    }
    my$ncount = () = $intronseq =~ m/N/g;
    if($ncount/$intronsize > 0.75){
      push(@comments, "$ncount Ns in fragment of size $intronsize: repeat interruption\n");
      $return = "1";
    }
  }
  close SEQFILE;
  return $return;
}

# split_gtf_for_optm splits large GTF files in files of max 5000 lines
# it keeps the genes intact (splits are only between genes) and
# returns the new filenames. 
# called from main

sub split_gtf_for_optm{

	my($inf)=(@_);
	my$sub="1";
	my$cutoff=5000;
	my@files=();
	my$splitflag;
	(my$filename=$inf) =~ s/\/.*\///;

# keep splitting the file until all files contain less than 5000 lines
	while(1){
		my$wc=`wc -l $inf`;
		(my$linenr)=$wc=~/(\d+)\s/;
		if($linenr<($cutoff+1)){
# quit if file has proper size
			unless($splitflag){
				push(@files, $inf);
				return (@files);
			}
# if remainder of file is below size cutoff: rename file
# create new temporary directories for all new files
			my$nr = rand;
			$nr =~ s/^0\.//;
			my$tmpdir="/$localdir/tmp.$nr";
			graceful_die ("Unique dir $tmpdir exists") if ( -d "$tmpdir");
			mkdir "$tmpdir" || graceful_die ("cannot create temporary directory $tmpdir: $!");
			(my$outf=$filename)=~ s/\.gtf$/\.$sub\.gtf/;
			system "mv $inf $tmpdir/$outf";
			push(@files, "$tmpdir/$outf");
			last;
		}else{
			$splitflag="1";
			my$tail=`tail +$cutoff $inf | head -n 1`;
			unless ($tail=~ /gene_id/){
				my$newcutoff=$cutoff+1;
				$tail =`tail +$newcutoff $inf | head -n 1`;
			}
			die "cannot find cutoff\n" unless ($tail =~ /gene_id/);
			(my$gene_id)=$tail =~ /gene_id \"(.*?)\"/;

# open file and read until gene_id
			my@newfile=();
			my@oldfile=();
			my$gotcha;
			open(GTF, "$inf") || die "cannot open $inf: $!\n";
			while (<GTF>){
				if(($_ =~ /gene_id \"$gene_id\"/) || ($gotcha)){
					push(@oldfile, $_);
					$gotcha="1";
				}else{
					push(@newfile, $_);
				}
			}
			close GTF;
# create new temporary directories for all new files
			my$nr = rand;
			$nr =~ s/^0\.//;
			my$tmpdir="/$localdir/tmp.$nr";
			graceful_die ("Unique dir $tmpdir exists") if ( -d "$tmpdir");
			mkdir "$tmpdir" || graceful_die ("cannot create temporary directory $tmpdir: $!");
			(my$outf=$filename)=~ s/\.gtf$/\.$sub\.gtf/;
			open(OUTF, ">$tmpdir/$outf");
			print OUTF @newfile;
			close OUTF;
			push(@files, "$tmpdir/$outf");
# put shorter version of file in place of original file
			open(OUTF2, ">$inf");
			print OUTF2 @oldfile;
			close OUTF2;
		}
		$sub++;
	}
	return(@files);
}

# remove_small removes all sequences <4 aa. Blast will not accept these as input
# called from run_synteny_method

sub remove_small{
  my($inputfilename)=(@_);
  my@kept=();
  open (INF, $inputfilename) or graceful_die ("cannot open inputfile $inputfilename: $!");
  my $fasta = new FAlite(*INF);
  while(my $entry = $fasta->nextEntry) {
    my$def = $entry->def;
    my$seq = $entry->seq;
    if(length $seq > 3){
      push(@kept, "$def\n", "$seq\n");
    }
  }
  close INF;
  open (OUTF, ">$inputfilename");
  print OUTF @kept;
  close OUTF;
  return;

}

# create_informant_blastdb gets the syntenic sequences for blast vs putative pseudogenes
# called from run_synteny_method

sub create_informant_blastdb{
  my$blastseqfile = shift @_;
  my@informant_hits = @_;

# some regions are selected twice. Keep in hash to remove the duplicates
  my%done=();
# keep track of chromosomes
  my%chrs=();

  foreach(@informant_hits){
    my@line = split ("\t", $_);
    my$musstart = $line[1];
    my$musend = $line[2];
    if($musstart > $musend){
      ($musstart, $musend) = ($musend, $musstart);
    }
    next if (exists($done{"$musstart".'TO'."$musend"}));
    $done{"$musstart".'TO'."$musend"}="1";
    my$chrnr = $line[0];
    $chrs{$chrnr} = "1";
    my$outfile="$tmpdir/$chrnr".'.informant.gtf';
    my$genename=$chrnr.'_'.$musstart;
    open(GTF, ">>$outfile") || graceful_die ("cannot create $outfile"); 
    print GTF 
      "$chrnr\tINFO\tCDS\t$musstart\t$musend\t\.\t+\t\.\tgene_id \"$genename\"\; transcript_id \"$genename\"\;\n";
    close GTF;
  }

# for every created gtf file, get the informant sequence
  foreach my$gtf(sort keys(%chrs)){
    my$file="$tmpdir/$gtf".'.informant.gtf';
    my$seq='chr'.$gtf.'.fa';
    system "/bio/bin/gtf2fa.pl -allcap -nogtf $file $informant_seqdir/$seq >> $blastseqfile";
  }

# create blastdb from outputfile (store STDERR info in unused array)
  my@comments = `xdformat -n -q3 $blastseqfile 2>&1`;
}

# format_params reads and validates the input from the parameter inputfile
# called from main

sub format_params{
	(my$inf)=(@_);
	my%hs=();
	my@faults=();
	open (INF, "$inf" ) || die "cannot open parameter file $inf: $!\n";
	while(<INF>){
		next if ($_ =~ /\#/);
		next unless ($_ =~ /=/);
		$_=~ s/\s//g;
		chomp $_;
		my($key, $value) = split('=', $_);
		$hs{$key} = $value;
	}
	close INF;

	$BLAST 			= $hs{wublast};
	$BLASTP			= $BLAST.'p';
	$BLASTN			= $BLAST.'n';
	($TBLASTN		= $BLAST)=~s/blast/tblastn/;
	$SIM4 			= $hs{sim4};
	$EST2GENO 		= $hs{est2genome};
	$parent_blastfile	= $hs{parent_blastfile};
	$parent_geno_seqs	= $hs{parent_geno_seqs};
	$parent_clone_seqs	= $hs{parent_clone_seqs};
	$target_seqdir		= $hs{target_seqdir};
	$informant_seqdir	= $hs{informant_seqdir};
	$prot_db		= $hs{prot_db};
	$prot_locations_file 	= $hs{prot_locations_file};
	$protein_sequences	= $hs{protein_sequences};
	$prot_loc_inputdir	= $hs{prot_loc_inputdir};
	$synteny_dir		= $hs{synteny_dir};
	$localdir		= $hs{tmpdir};
	$species		= $hs{species};
	$run_gene_pred		= $hs{run_gene_pred};

	# test existence of blast and alignment programs
	unless( -e $BLASTP ){
		push (@faults,  "cannot find $BLASTP in $BLAST");
	}
	unless( -e $BLASTN ){
		push (@faults,  "cannot find $BLASTN in $BLAST");
	}
	unless( -e $TBLASTN ){
		push (@faults,  "cannot find $TBLASTN in $BLAST");
	}

	$alignment_program="sim4";
	if( $SIM4 ){
		unless( -e $SIM4 ){
			push (@faults,  "cannot find $SIM4, using EST_GENOME instead");
			$alignment_program="est2geno";
		}
	}else{
		$alignment_program="est2geno";
	}

	unless( -e $EST2GENO ){
		push (@faults,  "cannot find $EST2GENO");
	}

	# see if all blast databases are available
	my@extensions=('xnd', 'xns', 'xnt');
	foreach(@extensions){
		my$findfile=$parent_blastfile.".$_";
		unless( -s $findfile ){
			push (@faults,  "cannot find file $findfile for $parent_blastfile");
		}
	}
	@extensions=('xpd', 'xps', 'xpt');
	foreach(@extensions){
		my$findfile=$prot_db.".$_";
		unless( -s $findfile ){
			push (@faults,  "cannot find file $findfile for $prot_db");
		}
	}

# check necessary files and directories for intron alignment method only if it's needed

	unless($opt_o){
# see if the directories exist
		unless ( -d $parent_geno_seqs ){
			push (@faults,  "cannot find $parent_geno_seqs");
		}
		unless ( -d $parent_clone_seqs ){
			push (@faults,  "cannot find $parent_clone_seqs");
		}
	}


	unless ( -d "$target_seqdir" ){
		push (@faults,  "cannot find target dir $target_seqdir");
	}
	unless ( -e "$informant_seqdir/chr1.fa" ){
		push (@faults,  "cannot find informant dir $informant_seqdir/chr1.fa");
	}
# check necessary files and directories for conserved_synteny method only if it's needed
	unless($opt_p){
		unless ( -d $prot_loc_inputdir ){
			push (@faults,  "cannot find $prot_loc_inputdir");
		}
		unless ( -d $synteny_dir ){
			push (@faults,  "cannot find $synteny_dir");
		}
		unless ( -d $protein_sequences ){
			push (@faults,  "cannot find $protein_sequences");
		}
# test if this directory is writable
		open ("OUTF" , ">$protein_sequences/testfile") || push (@faults,  "$protein_sequences is not writable!");
		close OUTF;
		unlink "$protein_sequences/testfile";
	}
	unless ( -d $localdir ){
		push (@faults,  "cannot find $localdir");
	}
	# test if this directory is writable
	open ("OUTF" , ">$localdir/testfile") || push (@faults,  "$localdir is not writable!");
	close OUTF;
	unlink "$localdir/testfile";

	my$twinscan_driver = $hs{twinscan_driver};
	if($hs{run_gene_pred} eq "twinscan" || $hs{run_gene_pred} eq "iscan"){
		my$genome_params = $hs{genome_params};
		if($hs{run_gene_pred} eq "iscan"){
			$genome_params = $hs{iscan_parameters};
		}else{
			unless ( -s $genome_params ){
				push (@faults,  "$genome_params doesn't exist or is empty");
			}
			my$cons_params = $hs{cons_params};
			unless ( -s $cons_params ){
				push (@faults,  "$cons_params doesn't exist or is empty");
			}
			$genome_params .= " $cons_params";
		}
		unless ( -s $twinscan_driver ){
			push (@faults,  "$twinscan_driver doesn't exist or is empty");
		}
		$informant_blastdb = $hs{informant_blastdb};
		@extensions=('xnd', 'xns', 'xnt');
		foreach(@extensions){
			my$findfile=$informant_blastdb.".$_";
			unless( -s $findfile ){
				push (@faults,  "cannot find file $findfile for $informant_blastdb");
			}
		}
# create twinscan_driver command

		$TWINSCAN_COMMAND = "$twinscan_driver $hs{extraparams} ". 
    			"-a $hs{blastparams} $genome_params $informant_blastdb "; 

	}elsif($hs{run_gene_pred} eq "nscan"){
		$nscan_params = $hs{nscan_param_file};
		unless ( -s $nscan_params ){
			push (@faults,  "$nscan_params doesn't exist or is empty");
		}
		$iscan_params = $hs{iscan_param_file};
		unless ( -s $iscan_params ){
			push (@faults,  "$iscan_params doesn't exist or is empty");
		}
		$align_ext = $hs{multispec_align_extension};
		$align_path = $hs{align_path};
		$TWINSCAN_COMMAND = "$twinscan_driver $nscan_params $iscan_params ";
	}
	return @faults;
}

