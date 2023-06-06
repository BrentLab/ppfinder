#!/usr/bin/perl -w
 
use strict;
use LWP::UserAgent;
use Getopt::Std;

use vars qw($opt_f);
getopts('f:');


my$usage="$0 <chr> <genome version> <informant genome version>

Program to download \'net synteny\' tables from ucsc and convert them to
tables with blocks of synteny (instead of long syntenic regions with gaps) 
genome version can be found on the browser homepage, eg hg16 for human and mm4 for mouse
Program converts 'top' and 'gap' lines in a ucsc file to blocks of orthologous regions.

Options:	-f <file> take inputfile instead of doing download

";

my($nr, $genome, $informant) = @ARGV;
die $usage unless $informant;

$nr =~ s/^chr//;

my@file = ();

if(defined($opt_f)){
  my$file = $opt_f;
  open(INF, $file) || die "cannot open inputfile: $!";
  @file = <INF>;
  close INF;
}else{

# for some weird reason, the informant must start with a capital
  $informant =~ s/^(.)/\U$1/;

  my $browser = LWP::UserAgent->new;
  my $response = $browser->post(
    'http://hgwdev.cse.ucsc.edu/cgi-bin/hgText',
#http://genome.ucsc.edu/cgi-bin/hgText
#    'http://genome.ucsc.edu/cgi-bin/hgText',
 [
  "hgsid" => "30333052",
  "db" => "$genome",
  "table" => "$genome.net$informant",
  "position" => "chr$nr",
  "tbPosOrKeys" => "pos",
  "outputType" => "Tab-separated\x2C Choose fields...",
  "cmp_bin" => "ignored",
  "pat_bin" => "",
  "cmp_level" => "ignored",
  "pat_level" => "",
  "dd_tName" => "does",
  "pat_tName" => "*",
  "cmp_tStart" => "ignored",
  "pat_tStart" => "",
  "cmp_tEnd" => "ignored",
  "pat_tEnd" => "",
  "dd_strand" => "does",
  "pat_strand" => "*",
  "dd_qName" => "does",
  "pat_qName" => "*",
  "cmp_qStart" => "ignored",
  "pat_qStart" => "",
  "cmp_qEnd" => "ignored",
  "pat_qEnd" => "",
  "cmp_chainId" => "ignored",
  "pat_chainId" => "",
  "cmp_ali" => "ignored",
  "pat_ali" => "",
  "cmp_score" => "ignored",
  "pat_score" => "",
  "cmp_qOver" => "ignored",
  "pat_qOver" => "",
  "cmp_qFar" => "ignored",
  "pat_qFar" => "",
  "cmp_qDup" => "ignored",
  "pat_qDup" => "",
  "dd_type" => "does",
  "pat_type" => "*",
  "cmp_tN" => "ignored",
  "pat_tN" => "",
  "cmp_qN" => "ignored",
  "pat_qN" => "",
  "cmp_tR" => "ignored",
  "pat_tR" => "",
  "cmp_qR" => "ignored",
  "pat_qR" => "",
  "cmp_tNewR" => "ignored",
  "pat_tNewR" => "",
  "cmp_qNewR" => "ignored",
  "pat_qNewR" => "",
  "cmp_tOldR" => "ignored",
  "pat_tOldR" => "",
  "cmp_qOldR" => "ignored",
  "pat_qOldR" => "",
  "cmp_tTrf" => "ignored",
  "pat_tTrf" => "",
  "cmp_qTrf" => "ignored",
  "pat_qTrf" => "",
  "log_rawQuery" => "AND",
  "rawQuery" => "",
  "origPhase" => "Get results",
  "boolshad.field_bin" => "1",
  "field_level" => "on",
  "boolshad.field_level" => "1",
  "field_tName" => "on",
  "boolshad.field_tName" => "1",
  "field_tStart" => "on",
  "boolshad.field_tStart" => "1",
  "field_tEnd" => "on",
  "boolshad.field_tEnd" => "1",
  "boolshad.field_strand" => "1",
  "field_qName" => "on",
  "boolshad.field_qName" => "1",
  "field_qStart" => "on",
  "boolshad.field_qStart" => "1",
  "field_qEnd" => "on",
  "boolshad.field_qEnd" => "1",
  "boolshad.field_chainId" => "1",
  "boolshad.field_ali" => "1",
  "boolshad.field_score" => "1",
  "boolshad.field_qOver" => "1",
  "boolshad.field_qFar" => "1",
  "boolshad.field_qDup" => "1",
  "field_type" => "on",
  "boolshad.field_type" => "1",
  "boolshad.field_tN" => "1",
  "boolshad.field_qN" => "1",
  "boolshad.field_tR" => "1",
  "boolshad.field_qR" => "1",
  "boolshad.field_tNewR" => "1",
  "boolshad.field_qNewR" => "1",
  "boolshad.field_tOldR" => "1",
  "boolshad.field_qOldR" => "1",
  "boolshad.field_tTrf" => "1",
  "boolshad.field_qTrf" => "1",
  "phase" => "Get these fields",
 ]
);

  die "Error: ", $response->status_line
   unless $response->is_success;

  my$content = $response->content;
  if($content =~ /No results/){
    print STDERR "No results for chr$nr\n";
    exit(0);
  }

  @file = split("\n", $content);
}

my%hs = ();
foreach my$line(@file){
  next if ($line =~ /^\#/);
  chomp $line;
  my@fields = split ("\t", $line);
  $hs{$fields[2]} = $line;
}
@file = ();
my@keys=sort numerically keys %hs;
foreach(@keys){
  push(@file, $hs{$_});
}


my($humchr_top, $humstart_top, $humend_top, $muschr_top, $musstart_top, $musend_top);
my($humchr_3, $humstart_3, $humend_3, $muschr_3, $musstart_3, $musend_3);
my($humchr_5, $humstart_5, $humend_5, $muschr_5, $musstart_5, $musend_5);
my($humchr_7, $humstart_7, $humend_7, $muschr_7, $musstart_7, $musend_7);
my($continue_top, $continue_3, $continue_5, $continue_7);
my$prevstart = "0";


foreach(@file){
  next if ($_ =~ /^\#/);
  chomp $_;

# top is always valid
# secondline is always gap
# if third level is NonSyn, the 4th level is gap for the 3rd level
# in fact, nonSyn is always 3, 5, or 7
# all gaps are in even numbers (max 8)
# inv is in level 3 or 5
# top is in level 1
# syn is a lineup on the same chromosome as the gap in the level above it. - dodgy at best?
# so we're looking at nested data. We could sort, but that could get messy
# maybe, whenever you see a top (level 1, 3, 5, 7), decide if it's big enough, then continue

  my@line = split("\t", $_);
  if($line[2] < $prevstart){
    die "unsorted file!";
  }
  $prevstart = $line[2];

  if($line[7] eq "top"){
    $continue_top = "0";

# if there is a remainder of the previous region, print it
    if($humend_top){
      print "$humchr_top\t$humstart_top\t$humend_top\t$muschr_top\t$musstart_top\t$musend_top\n";
      ($humchr_top, $humstart_top, $humend_top, $muschr_top, $musstart_top, $musend_top) = ("", "", "", "", "", "");
    }
# only continue if the region is big enough
    if($line[3] - $line[2] > 10000){
#    if($line[3] - $line[2] > 4500){
      $humchr_top = $line[1];
      $humstart_top = $line[2];
      $humend_top = $line[3];
      $muschr_top = $line[4];
      $musstart_top = $line[5];
      $musend_top = $line[6];
      $continue_top = "1";
      die "Start and end in wrong order" if ($humstart_top > $humend_top);
    }
  }elsif($line[0] eq "2"){
    if($continue_top){
      ($humstart_top, $musstart_top) = create_block($humstart_top, $musstart_top, $humchr_top, $muschr_top,  @line);
    }
  }elsif($line[0] eq "3"){
    $continue_3 = "0";

# if there is a remainder of the previous region, print it
    if($humend_3){
      print "$humchr_3\t$humstart_3\t$humend_3\t$muschr_3\t$musstart_3\t$musend_3\n";
      ($humchr_3, $humstart_3, $humend_3, $muschr_3, $musstart_3, $musend_3) = ("", "", "", "", "", "");
    }
# only continue if the region is big enough
    if($line[3] - $line[2] > 10000){
#    if($line[3] - $line[2] > 4500){
      $humchr_3 = $line[1];
      $humstart_3 = $line[2];
      $humend_3 = $line[3];
      $muschr_3 = $line[4];
      $musstart_3 = $line[5];
      $musend_3 = $line[6];
      $continue_3 = "1";
      die "Start and end in wrong order" if ($humstart_3 > $humend_3);
    }
  }elsif($line[0] eq "4"){
# only continue if the region is large enough
    if($continue_3){
      ($humstart_3, $musstart_3) = create_block($humstart_3, $musstart_3, $humchr_3, $muschr_3, @line);
    }
  }elsif($line[0] eq "5"){
    $continue_5 = "0";
# print previous line and reset
    if($humend_5){
      print "$humchr_5\t$humstart_5\t$humend_5\t$muschr_5\t$musstart_5\t$musend_5\n";
      ($humchr_5, $humstart_5, $humend_5, $muschr_5, $musstart_5, $musend_5) = ("", "", "", "", "", "");
    }
    if($line[3] - $line[2] > 10000){
#    if($line[3] - $line[2] > 4500){
      $humchr_5 = $line[1];
      $humstart_5 = $line[2];
      $humend_5 = $line[3];
      $muschr_5 = $line[4];
      $musstart_5 = $line[5];
      $musend_5 = $line[6];
      die "Start and end in wrong order" if ($humstart_5 > $humend_5);
      $continue_5 = "1";
    }
  }elsif($line[0] eq "6"){
    if($continue_5){
      ($humstart_5, $musstart_5) = create_block($humstart_5, $musstart_5, $humchr_5, $muschr_5, @line);
    }
  }elsif($line[0] eq "7"){
    $continue_7 = "0";
    if($humend_7){
      print "$humchr_7\t$humstart_7\t$humend_7\t$muschr_7\t$musstart_7\t$musend_7\n";
      ($humchr_7, $humstart_7, $humend_7, $muschr_7, $musstart_7, $musend_7) = ("", "", "", "", "", "");
    }
    if($line[3] - $line[2] > 10000){
#    if($line[3] - $line[2] > 45000){
      $humchr_7 = $line[1];
      $humstart_7 = $line[2];
      $humend_7 = $line[3];
      $muschr_7 = $line[4];
      $musstart_7 = $line[5];
      $musend_7 = $line[6];
      $continue_7 = "1";
      die "Start and end in wrong order" if ($humstart_7 > $humend_7);
    }
  }elsif($line[0] eq "8"){
    if($continue_7){
      ($humstart_7, $musstart_7) = create_block($humstart_7, $musstart_7, $humchr_7, $muschr_7, @line);
    }
  }else{
#    my$diff = $line[3] - $line [2];
#    print STDERR "do not recognize line starting with $line[0]\t$diff\n";
  }

}
# print all remaining regions

if($humend_top){
  print "$humchr_top\t$humstart_top\t$humend_top\t$muschr_top\t$musstart_top\t$musend_top\n";
}
if($humend_3){
  print "$humchr_3\t$humstart_3\t$humend_3\t$muschr_3\t$musstart_3\t$musend_3\n";
}
if($humend_5){
  print "$humchr_5\t$humstart_5\t$humend_5\t$muschr_5\t$musstart_5\t$musend_5\n";
}
if($humend_7){
  print "$humchr_7\t$humstart_7\t$humend_7\t$muschr_7\t$musstart_7\t$musend_7\n";
}


############ SUBROUTINES ###############

sub numerically {$a <=> $b}

sub create_block{
  my($humstart, $musstart, $humchr, $muschr,  @line) = (@_);

  die "level $line[0] is not a gap:$_" unless ($line[7] eq "gap");
  die "no comparison start found for @line" unless ($humstart =~ /\d/);
  die "Something wrong with $_: gap $muschr not part of $line[4]" unless ($line[4] eq $muschr);
  die "Start and end in wrong order" if ($line[2] > $line[3]);
  
# create new block
  my$newhumstart = $humstart;
  my$newhumend = $line[2]-1;	#one before start of gap        
  my$newmusstart = $musstart;
  my$newmusend = $line[5]-1;	#one before start of gap
  $musstart = $line[6]+1;	# new start is one after end of gap

  if($line[5] == $line[6]){		# start and end the same (it's unclear to which block the single nt belongs)
    $newmusend = $line[5];
    $musstart = $line[6];
  }
  if($line[2] == $line[3]){
    print STDERR "gap in human: $line[2]";
  }

  my$humdif = $newhumend - $newhumstart;
  my$musdif = $newmusend - $newmusstart;

  print "$humchr\t$newhumstart\t$newhumend\t$muschr\t$newmusstart\t$newmusend\n";
  $humstart = $line[3]+1;	#one after end of gap

  return ($humstart, $musstart); 

}

