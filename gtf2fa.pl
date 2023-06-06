#!/usr/bin/perl

# this program was created by Michael Stevens in the Brent Lab

use Getopt::Long;
use GTF;
use Data::Dumper;
use strict;

my %c = (
# Global
  norc      => '',    # don't reverse complement - strands
  tx_id     => '',    # Only process records matching this tx_id
  gid       => '',    # output gene id instead of transcript id
  h         => '',    # usage
# Parser
  exon      => '',    # parse exons not features (ie, ignore lines w/o exon in them. features=UTR5,start_codon,CDS,stop_codon,UTR3)
  cds       => '',    # only use CDS features (ie, ignore exon, UTR5,start_codon,stop_codon,UTR3)
  noGTF     => '',    # don't use gtf.pm, use inline parser
  nocheck   => '',    # disable GTF.pm's checking
# DNA mode
  tx        => '',    # print just the tx, ie verses the whole genomic seq
  nocap     => '',    #  don't uc CDS/exons (note: stop_exon isn't part of CDS) 
  allcap    => '',    # capitalize everything
# DNA's tx mode
  flank     => '',    # flank 1st and last exon by x bases, will report in header if flank runs off either end of fa 
# protein mode 
  ptx       => '',    # print out the protein sequence instead of DNA (not defined w/ -exon)
  split     => '',    # one fa entry per line in gtf file only works w/ ptx
  frame     =>'1',     # Use the frame in each line to determine frame, w/o this frame is calculated from the start codon, always use frame in cds
);

my $usage = "
$0 [options] <gtf file> <fa file>

Given pull out the genomic seq from a fa file specified by a gtf file.  Assumes just one fa entry per fasta file.  
By default, it ignores lines with exon in them.

Options:
@{ [ qx|perl -ne 'next unless (\$s or m/\%c/);last if m/^\\)/;print \$_ if \$s++' $0|  ] }
";


$Data::Dumper::Varname = "Conf"; 
GetOptions(\%c,'norc','exon','cds','split','ftr=s','tx','ptx', 'nocap','allcap','gid',"flank=s",'tx_id=s','noGTF','nocheck','frame', "h") 
  || die join "\n",$usage, Dumper(\%c),'';

die join "\n", $usage, Dumper(\%c) if ($c{h} || @ARGV < 1);

#### validation ####
$c{cds} = 1 if ($c{ptx}); # 050906 to make readout_ptx works right when UTR are present

###################


my ($gtf_fn,$fa_fn) = @ARGV;

#
my $gtf = gtf_sort(gtf_hash($gtf_fn));

get_seq($gtf,$fa_fn);

sub comma { my $v = shift; 1 while ($v =~ s/([-+]?\d+)(\d{3})/$1,$2/); return $v; }

# assumes just one entry in fa file
sub get_seq {
  my ($gtf,$fn) = @_;
  open F, $fn or die "can't open $fn: $!";
  # idx is 1-based
  my ($pos,$idx);
  #warn "gtf has:",scalar(keys %$gtf)," tx's";
  my $ftr = $c{exon} ? 'exn' : 'ftr';
  # sort xscripts by 1)1st ftr bgn and 2) 1st ftr end
  for my $tx (
    sort {$gtf->{$a}{$ftr}[0][0] <=> $gtf->{$b}{$ftr}[0][0] 
      || $gtf->{$a}{$ftr}[-1][1] <=> $gtf->{$b}{$ftr}[-1][1]
    } keys %$gtf
  ) {

    # $ftr for this could be empty
    next unless (exists($gtf->{$tx}{$ftr}) && @{$gtf->{$tx}{$ftr}[0]});

    # flank is 0 unless specified
    my $tx_start = $gtf->{$tx}{$ftr}[0][0] - $c{flank};
    # in case flank runs of start
    $tx_start = 1 if $tx_start < 0;
    while (<F>) {
      next if (/^#/ || /^\s*$/ || /^>/);
      # idx is last char read so far
      if ($idx+length($_)-1 >= $tx_start ) {
        # warn "hitme $tx @ ",comma($gtf->{$tx}{$ftr}[0][0]);
        # set to read from idx (+1 for newline
        seek(F, $tx_start - ($idx+length($_)+1),1);
        my $pos = tell(F);
        ($c{ptx}) ? readout_ptx(*F,$tx,$gtf->{$tx}) : readout_tx(*F,$tx,$gtf->{$tx});
        # reset to idx to continue with next tx ( may have same start)
        seek(F,$pos,0);
        $idx = $tx_start-1;
        last;  #go to next gtf
      } else {
        $idx += (length($_)-1);
      }

    }
  }
  close F;
}


# 050906
# doesn't work correctly if there's UTR in the ftr
# instead of fixing it, just ensure -cds if set if -ptx is set
# which should avoid the problem, maybe fix this later
sub readout_ptx {
  my ($fh,$tx_id,$tx) = @_;

  my %AA = ( 
    ATT=>'I', ATC=>'I', ATA=>'I',
    CTT=>'L', CTC=>'L', CTA=>'L', CTG=>'L', TTA=>'L', TTG=>'L',
    GTT=>'V', GTC=>'V', GTA=>'V', GTG=>'V',
    TTT=>'F', TTC=>'F',
    ATG=>'M',
    TGT=>'C', TGC=>'C',
    GCT=>'A', GCC=>'A', GCA=>'A', GCG=>'A',
    GGT=>'G', GGC=>'G', GGA=>'G', GGG=>'G',
    CCT=>'P', CCC=>'P', CCA=>'P', CCG=>'P',
    ACT=>'T', ACC=>'T', ACA=>'T', ACG=>'T',
    TCT=>'S', TCC=>'S', TCA=>'S', TCG=>'S', AGT=>'S', AGC=>'S',
    TAT=>'Y', TAC=>'Y',
    TGG=>'W',
    CAA=>'Q', CAG=>'Q',
    AAT=>'N', AAC=>'N',
    CAT=>'H', CAC=>'H',
    GAA=>'E', GAG=>'E',
    GAT=>'D', GAC=>'D',
    AAA=>'K', AAG=>'K',
    CGT=>'R', CGC=>'R', CGA=>'R', CGG=>'R', AGA=>'R', AGG=>'R',
    #Stop codons     Stop    TAA, TAG, TGA
    TAA=>'*', TAG=>'*', TGA=>'*',
  );


  my $strand = $tx->{strand};

  # ptx _always_ use CDS
  my $ftr = $tx->{ftr};
  # 1st check whether gtf contains CDS features
  return unless (grep { $_->[2] eq 'CDS' } @$ftr); 

  my $start = (grep { $_->[2] eq 'CDS' } @$ftr)[0]->[0]; 
  my $stop  = (grep { $_->[2] eq 'CDS' } @$ftr)[-1]->[1];
  my $len;

#  my ($stop) = , $stop, $len) = ($ftr->[0][0],$ftr->[-1][1],0);
#  my ($start, $stop, $len) = ($tx->{start_codon},$tx->{stop_codon}+3, 0);


  my $tlen = 1+ $stop - $start; 

#printf "readout: >$tx_id 1+%s-%s=%s\n",$stop,$start,$tlen;
  my $id = $c{gid} ? $tx->{gid} : $tx_id; 
  my $hdr = ">$id start:$start end:$stop strand $strand\n";

# read genomic seq
  my $seq;
  while (<$fh>) {
    next if (/^#/ || /^\s*$/ || /^>/);
    chomp;  
    my $f = $len + length($_);
    if ($f < $tlen) {
      $seq .=  $_; 
      $len  = $f;
    } else {
      #   print "\n$f,$len,$tlen",$tlen-$len,"\n";
      $seq .= substr($_,0, ($tlen - $len)); 
      last;   
    }
  }

  my $base = $ftr->[0][0];
  my $tx_fa;

  my @splice = (0);
  my @coord;
  my @frame;
  for (@$ftr) {
    # start_codons are redundant 
    next unless $_->[2] eq 'CDS';
    push @frame, $_->[3];

    # susbstr is 0-based, gtf is 1-based
    $tx_fa .= substr($seq,$_->[0] - $base ,1+$_->[1]-$_->[0]);
    push @splice, length($tx_fa);
    push @coord, [$_->[0],$_->[1]];
  }

  $seq = uc($tx_fa);

  # rc it?
  if (!$c{norc} && $strand eq '-') {
    $seq =~ tr/atgcATGC/tacgTACG/;
    $seq = reverse(scalar($seq));
    @coord = reverse(@coord);
    @frame = reverse(@frame);
    @splice = map { length($tx_fa) - $_ } reverse(@splice);
  }


  if ( $c{split} ) {
    my $frame = 0;
   for (1..$#splice) {

     # explicitly get frame from line  
     $frame =  $frame[$_-1] if ($c{frame});
       
       
     print join "\t", ">$id.exon_".$_,
     #                     "source=",
                      "source_id=$id",
                      "start=".$coord[$_-1][0],
                      "end=".$coord[$_-1][1] ;
     print "\n";
     my $subseq = substr($seq, $splice[$_-1], $splice[$_]-$splice[$_-1]); 
#     $subseq = substr($subseq,3-$frame) if ($frame);
#     $frame = length($subseq)%3;
#     $subseq =substr($subseq,0,-$frame) if ($frame);
      frame_trim(\$subseq,$frame);
     $subseq =~ s/(...)/$AA{$1}/g;
     print $subseq,"\n";
   }
  } else {
   frame_trim(\$seq, $frame[0]);
   $seq =~ s/(...)/$AA{$1}/g;
  print $hdr,$seq,"\n" if ($seq);
  }

}

sub frame_trim {
my ($seq, $frame) = @_;
     $$seq = substr($$seq,$frame);
     $frame = length($$seq)%3;
     $$seq =substr($$seq,0,-$frame) if ($frame);
}

# features are sorted by start
sub readout_tx {
  my ($fh,$tx_id,$tx) = @_;
  
  my $strand = $tx->{strand};
  my $ftr = $c{exon} ? $tx->{exn} : $tx->{ftr};
  my ($start, $stop, $len) = ($ftr->[0][0],$ftr->[-1][1],0);
  my ($lflank,$rflank) = ($c{flank},$c{flank});

  if ($start > $lflank) {
    # flank is 0 unless given
    $start -= $lflank;
  } else { # not enough file
     $lflank = $start -1;
     $start = 1;
  }

  $stop  += $rflank;
  my $tlen = 1+ $stop - $start;

#printf "readout: >$tx_id 1+%s-%s=%s\n",$stop,$start,$tlen;
  my $id = $c{gid} ? $tx->{gid} : $tx_id;

#  # if c{exon} then the previous setting are correct
#  if ($c{tx} && !$c{exon}) {
#    # assumes ftr is sorted, then grab the 1st and last CDS
#    $start = ( $tx->{start_codon} || (grep { $_->[2] eq 'CDS' } @$ftr)[0]->[0] ); 
#    $stop  = ( $tx->{stop_codon} || (grep { $_->[2] eq 'CDS' } @$ftr)[-1]->[1] ); 
#
#  }

# read genomic seq
  my $seq;
  while (<$fh>) {
    next if (/^#/ || /^\s*$/ || /^>/);
    chomp;
    my $f = $len + length($_);
    if ($f < $tlen) {
      $seq .=  $_;
      $len  = $f;
    } else {
      #   print "\n$f,$len,$tlen",$tlen-$len,"\n";
      $seq .= substr($_,0, ($tlen - $len));
      last;
    }
  }

  # ran out of file bfr reaching tlen 
  #(if not rflank, then gtf shouldn't be going off end, should add to validation, not here)
  if ($rflank && length($seq) < $tlen ) { $stop = $start + length($seq) -1;$rflank = $stop - $ftr->[-1][1]   }

  my $hdr = ">$id start:$start end:$stop strand $strand";
  # only report flanks if clipped
  $hdr .= " (lflank = $lflank) " if ($lflank < $c{flank});
  $hdr .= " (rflank = $rflank) " if ($rflank < $c{clank});
  $hdr .= "\n";

  # print genomic or tx
  if ($c{tx} || !$c{nocap}) {

    #1st lower it
    $seq = lc($seq); # unless ($c{tx});

    my $tx_fa = substr($seq,0,$lflank);
    my $base = $ftr->[0][0] - $lflank;

    for (@$ftr) {
      # start_codons are redundant 
      next if $_->[2] eq 'start_codon';

      # susbstr is 0-based, gtf is 1-based
      my $str = substr($seq,$_->[0] - $base ,1+$_->[1]-$_->[0]);

      # stop codon, not seq but include in lc for when using flank
      #if (!$c{nocap} && $_->[2] ne 'stop_codon') {
     # Nope: current consensus is to cap stop_codon
     if (!$c{nocap}) {
       $str = uc($str); # unless ($c{nocap});
        # if genomic put uc ftr back
        substr($seq,$_->[0] - $base ,1+$_->[1]-$_->[0]) = $str unless ($c{tx});
      }

      $tx_fa .= $str if ($c{tx});
    }

    if ($c{tx}) {
      $seq = $rflank ? $tx_fa.substr($seq,-1*$rflank) : $tx_fa;
    }
  } 

  # rc it?
  if (!$c{norc} && $strand eq '-') {
    $seq =~ tr/atgcATGC/tacgTACG/;
    $seq = reverse(scalar($seq));
  }

  # forget all the work we've done, just capitalize everything
  $seq = uc($seq) if ($c{allcap}); 
  print $hdr,$seq,"\n" if ($seq);

}

sub gtf_hash {
  my ($fn) = @_;
  my %rtn;

  if ($c{noGTF}) {
    open F, $fn or die "can't open $fn: $!";
    while (<F>) {
      next if (/^#/ || /^\s*$/);
      chomp;
      my @line = split /\t/;
      my $tx_id = ($line[8] =~ /transcript_id "([^"]+)"/)[0];
      $rtn{$tx_id}{gid} = ($line[8] =~ /gene_id "([^"]+)"/)[0];
      $rtn{$tx_id}{chr}         = $line[0];
      $rtn{$tx_id}{strand}      = $line[6];
      $rtn{$tx_id}{start_codon} = $line[3] if $line[2] =~ /start_codon/;
      $rtn{$tx_id}{stop_codon}  = $line[3] if $line[2] =~ /stop_codon/;

      if ($line[2] =~ /exon/i) {
        push @{$rtn{$tx_id}{exn}}, [$line[3],$line[4],$line[2], $line[7]];
      } else {
       push @{$rtn{$tx_id}{ftr}}, [$line[3],$line[4],$line[2], $line[7]] if ($line[2] =~ /cds/i || !$c{cds}); 
      }

      #      ($chr,$start,$transid) = (split /\t/,$_)[0,3,8];
    }
    close F;
  } else { # use GTF

    my $tscripts = (GTF::new({ gtf_filename=> $fn, no_check=>$c{nocheck} }))->transcripts;
    for my $tx (@$tscripts) {
      my %gtf;
      $gtf{gid}    = $tx->gene->id;
      $gtf{chr}    = $tx->gene->seqname;
      $gtf{strand} = $tx->gene->strand;
      $gtf{start_codon} = $tx->start_codons->[0]->{Start} if ($tx->start_codons->[0]);
      $gtf{stop_codon} = $tx->stop_codons->[0]->{Start} if ($tx->stop_codons->[0]);

      for my $cds (@{$tx->cds}) {
        push @{$gtf{ftr}}, [0+$cds->start,0+$cds->stop,$cds->type, $cds->frame];
      } 

      unless ($c{cds}) {
        for my $cds (@{$tx->utr5}, @{$tx->start_codons},@{$tx->stop_codons}, @{$tx->utr3}) {
          push @{$gtf{ftr}}, [0+$cds->start,0+$cds->stop,$cds->type, $cds->frame];
        } 
      }

      for my $xon (@{$tx->exons}) {
        push @{$gtf{exn}}, [0+$xon->start,0+$xon->stop, $xon->type, $xon->frame];
      }
      $rtn{$tx->id} = \%gtf;
    }       

  }
  #print Dumper(\%rtn);exit;

  %rtn = ( $c{tx_id}=>$rtn{$c{tx_id}} ) if ($c{tx_id});
  return \%rtn; 
}

sub gtf_sort {
  my ($gtf) = @_;
  # sort features in each tx by: 1) ftr start value 2) ftr stop val 
  for my $g (values %$gtf) {
    for ($g->{exn},$g->{ftr}){ 
      next unless $_;
      $_ = [sort {$a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] } @$_];
    }
  }
  # sort tx's by: 1) 1st exn start 2) last exn stop
  # $gtf = {
  #    sort {$a->{exn}[0][0] <=> $b->{exn}[0][0] || $a->{exn}[-1][1] <=> $b->{exn}[-1][1]} keys %$gtf
  #  }; 
  return $gtf;
}

#chr2    iscan   stop_codon      12902657        12902659        0.      +       0       gene_id "13.seq.masked.003"; transcript_id "13.seq.masked.003.1";
#chr2    iscan   CDS     238759069       238759191       0.      +       0       gene_id "239.seq.masked.004"; transcript_id "239.seq.masked.004.1";
#chr2    iscan   CDS     238761908       238761966       0.      +       0       gene_id "239.seq.masked.004"; transcript_id "239.seq.masked.004.1";
#chr2    iscan   start_codon     238761967       238761969       0.      +       0       gene_id "239.seq.masked.004"; transcript_id "239.seq.masked.004.1";
#chr1    A-list  exon    1322311 1322603 .       +       .       gene_id "NM_199121"; transcript_id "NM_199121"; exon_number "1"; exon_id "NM_199121.1";
#chr1    A-list  UTR     1322311 1322603 .       +       .       gene_id "NM_199121"; transcript_id "NM_199121"; exon_number "1"; exon_id "NM_199121.1";
#chr1    A-list  exon    1325863 1327547 .       +       .       gene_id "NM_199121"; transcript_id "NM_199121"; exon_number "3"; exon_id "NM_199121.3";
#chr1    A-list  UTR     1325863 1325867 .       +       .       gene_id "NM_199121"; transcript_id "NM_199121"; exon_number "3"; exon_id "NM_199121.3";
#chr1    A-list  CDS     1325868 1326569 .       +       .       gene_id "NM_199121"; transcript_id "NM_199121"; exon_number "3"; exon_id "NM_199121.3";
#chr1    A-list  UTR     1326570 1327547 .       +       .       gene_id "NM_199121"; transcript_id "NM_199121"; exon_number "3"; exon_id "NM_199121.3";
#chr1    A-list  start_codon     1325868 1325870 0       +       .       gene_id "NM_199121"; transcript_id "NM_199121"; exon_number "3"; exon_id "NM_199121.3";
#chr1    A-list  stop_codon      1326567 1326569 0       +       .       gene_id "NM_199121"; transcript_id "NM_199121"; exon_number "3"; exon_id "NM_199121.3";


# NOTES:
# 050309  added readout_ptx and ptx to print out protein tx
# 050309  added tx_id to choose one tx_id from the gtf
# 050309  added split to printout each exon on seperate line, only defn'd for ptx
# 050314  removed using start/stop_codon to set start and stop, old way seemed to the right, ie ftr->[0][0] - ftr->[-1][1]
# 050906 to make readout_ptx works right when UTR are present, always set -cds when -ptx
# 050906 make -frame=1 by default
