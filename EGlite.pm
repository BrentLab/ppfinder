package EGlite;
use strict;
use overload '""' => '_overload';
###############################################################################
# EGlite
###############################################################################
sub new {
	my ($class, $fh) = @_;
	if (ref $fh !~ /GLOB/)
		{die "EGlite error: new expects a GLOB reference not $fh\n"}
	my $this = bless {};
	$this->{FH} = $fh;
	$this->{LASTLINE} = "";
	if ($this->_parseHeader) {
		# there are alignments
		$this->{REPORT_DONE} = 0;
		die if not defined $this->{GENOME};
		die if not defined $this->{EST};
		if ($this->_parseAlignment) {
			# there is alignment
			$this->{REPORT_DONE} = 1;
		} else {
			return 0;
		}
	}
	else {
		# empty report
		$this->{REPORT_DONE} = 1; 
		return 0;
	}
	return $this;
    
}
sub genome          {shift->{GENOME}}
sub est             {shift->{EST}}
sub score           {shift->{SCORE}}
sub pi              {shift->{PI}}
sub alignment       {shift->{ALIGNMENT}}
sub gb              {shift->{GENOME_BEGIN}}
sub ge              {shift->{GENOME_END}}
sub eb              {shift->{EST_BEGIN}}
sub ee              {shift->{EST_END}}
sub genomeBegin     {shift->{GENOME_BEGIN}}
sub genomeEnd       {shift->{GENOME_END}}
sub estBegin        {shift->{EST_BEGIN}}
sub estEnd          {shift->{EST_END}}
sub genomeDirection {shift->{GEN_DIRECTION}}
sub estDirection    {shift->{EST_DIRECTION}}
sub geneOrientation {shift->{SPLICE_DIRECTION}}

sub countElements   {
	my ($this) = @_;
	my $elements = $this->{ELEMENTS};
	my @elements = @$elements;
	return $#elements + 1;
}

sub hasMoreElements {
	my ($this) = @_;
	my $buf = $this->{ELEMENTSBUFFER};
	my @buf = @$buf;
	if( defined $buf[0]){
		return 1;
	}else{
		return 0;
	}
}

sub hasIntron {
	my ($this) = @_;
	my $elements = $this->{ELEMENTS};
	if ($#$elements <= 1) {
		return 0;
	} else {
		return 1;
	}
}

sub nextElement {
	my ($this) = @_;
#	$this->_fastForward or return 0;
	my $buffer = $this->{ELEMENTSBUFFER};	
	my $element = shift (@$buffer);
	return 0 unless defined $element;
	$this->{ELEMENTSBUFFER} = $buffer;
	return $element;
}

sub resetElements {
	my ($this) = @_;
	my $elements = $this->{ELEMENTS};
	my @elements = @$elements;
	my @buffer   = @elements;
	$this->{ELEMENTSBUFFER} = \@buffer;
#print $this->{ELEMENTSBUFFER}, "\n", $this->{ELEMENTS}, "\n";
	return 0 unless defined \@buffer;
	return 1;
}

sub _parseHeader {
	my ($this) = @_;
	my $FH = $this->{FH};
	my @introns = ();
	my @exons = ();
	my @elements = ();
	my @buffer = ();

	while(<$FH>) {
		if ($_ =~ /Note Best alignment is between (reversed|forward) est and (reversed|forward) genome, (but|and) splice\s+sites imply\s+(forward gene|REVERSED GENE)/ ) {
			$this->{EST_DIRECTION} = $1;
			$this->{GEN_DIRECTION} = $2;
			if ($4 eq "forward gene") {
				$this->{SPLICE_DIRECTION} = '+';
			} elsif ($4 eq "REVERSED GENE") {
				$this->{SPLICE_DIRECTION} = '-';
			} else {
				die "Unknown gene orientation found in alignment: $4\n";
			}
				
		}
		elsif ($_ =~ /^Exon\s/)    {
			my (undef, $score, $pi, $gb, $ge, $genome, $eb, $ee, $est) = split (/\s+/);
			my $exon = EGlite::Exon::new( GENOME=>$genome, 
							EST=>$est, 
							GB=>$gb, 
							'GE'=>$ge, 
							EB=>$eb, 
							EE=>$ee, 
							SCORE=>$score, 
							PI=>$pi, 
							DIRECTION=>".",
							PARENT=>$this);
			push @exons, $exon;
			push @elements, $exon;
		}
		elsif ($_ =~ /^Span\s/)    {
			my (undef, $score, $pi, $gb, $ge, undef, $eb, $ee, undef) = split (/\s+/);
			$this->{GENOME_BEGIN} = $gb;
			$this->{GENOME_END}   = $ge;
			$this->{EST_BEGIN}    = $eb;
			$this->{EST_END}      = $ee;
			$this->{SCORE}        = $score;
			$this->{PI}           = $pi;
		}
		elsif ($_ =~ /^(.)Intron\s/)    {
			my $direction = $1;
			my (undef, $score, $pi, $gb, $ge, $genome) = split (/\s+/);
			my $intron = EGlite::Intron::new( GENOME=>$genome, 
							GB=>$gb, 
							'GE'=>$ge, 
							SCORE=>$score, 
							PI=>$pi, 
							DIRECTION=>$direction,
							PARENT=>$this);
			push @introns, $intron;
			push @elements, $intron;
		}
		elsif ($_ =~ /^(.+)\s+vs\s+(.+):/) {
			$this->{LASTLINE} = $_; 
			$this->{GENOME}   = $1;
			$this->{EST}      = $2;
			$this->{EXONS} = \@exons;
			$this->{INTRONS} = \@introns;
			$this->{ELEMENTS} = \@elements;
			$this->resetElements; # Sets ELEMENTBUFFER
			return 1;
		}
		elsif ($_ =~ /^Alignment Score: (\d+)/) {
			$this->{LASTLINE} = $_;
			$this->{REPORT_DONE} = 1;
			$this->{SCORE} = $1;
			return 0; # there's nothing in the report
		}
	}
	return 0 unless (defined $this->{GENOME} && defined $this->{EST});
}

sub _parseAlignment {
	my ($this) = @_;
	return 0 if $this->{ALIGN_ALL_PARSED};

	my $est = $this->{EST};
	my $genome = $this->{GENOME};
	my $FH  = $this->{FH};
	
	#######################
	# get alignment lines #
	#######################
	
	my @algnline;
	while(<$FH>) {
		die if /^BLAST/;
#print "$_";
		if ($_ =~ /^WARNING:|^NOTE:|^ERROR:|^FATAL:/) {
			while(<$FH>) {last if $_ !~ /\S/}
		}
		elsif ($_ !~ /\S/)            {next}
		elsif ($_ =~ /^EXIT CODE/) {last}
		elsif ($_ =~ /^Alignment\s+Score:\s+([0-9\-\.]+)\s+$/)   {
			$this->{LASTLINE} = $_;
			$this->{SCORE}    = $1;
			$this->{ALIGN_ALL_PARSED} = 1;
#print "In Last Line\n";
			last;
		}
		#elsif ($_ =~ /^\s*($genome)\s+\d+\s+/) {
		elsif ($_ =~ /^\s*(\S+)\s+\d+\s+/) {
#print "In genome";
		    if ($1 eq $genome){
			push @algnline, $_;           # store the genome line
			my $l1 = <$FH>;              # either alignment line or est line
			#if ($l1 =~ /^\s*($est)\s+\d+\s+/) {
			 if ($l1 =~ /^\s*(\S+)\s+\d+\s+/) {   
			     if ($1 eq $est){
				push @algnline, "";  # dummy line, this is a -noseq option
				push @algnline, $l1; # so store a fake alignment and real sbjct
				next;
			    }
			}
			push @algnline, $l1;                 # grab/store the alignment line
			my $l2 = <$FH>; push @algnline, $l2; # grab/store the est line
		    }
		}
	}
	
	#########################
	# parse alignment lines #
	#########################
	my ($gl, $el, $as) = ("", "", "");
	my ($gb, $ge, $eb, $ee) = (0,0,0,0);
	my (@GL, @EL, @AS); # for better memory management
			
	for(my $i=0;$i<@algnline;$i+=3) {
#print "In..";
		#warn $algnline[$i], $algnline[$i+2];
		#$algnline[$i]   =~ /^\s*($genome)\s+(\d+)\s+(\S+)\s+(\d+)/;
	        $algnline[$i]   =~ /^\s*(\S+)\s+(\d+)\s+(\S+)\s+(\d+)/;
		if ($1 ne $genome){ return 0;}
		return 0 unless defined $4;
		$gl = $3; $gb = $2 unless $gb; $ge = $4;
		my $number = $2; # Catch the est_genome glitch that has 6 char width for coordinates
		
		# Found a bug where there was one base "C" in the alignment
		# and the name of the sequence had a C in it too. So index()
		# picked up the C in the name. So I have to make sure that 
		# index() starts after $gb

		my $offset = index($algnline[$i], $gl, index($algnline[$i], $number, 0));

		# est_genome has hardcoded 6 char width for coordinates which breaks here.
		if (length($number) > 6) {
			$offset -= (length($number) - 6);
		}
		$as = substr($algnline[$i+1], $offset, CORE::length($gl))
			if $algnline[$i+1];
		
		$algnline[$i+2] =~ /^\s*(\S+)\s+(\d+)\s+(\S+)\s+(\d+)/;
		if ($1 ne $est){ return 0;}
		return 0 unless defined $4;
		$el = $3; $eb = $2 unless $eb; $ee = $4;
		
$gl =~ s/[^ATGCNatgcn\-\.]//;
		push @GL, $gl; push @EL, $el; push @AS, $as;
#print "PARSER:\n$gl\n$as\n$el\n";
	}

	########################
	# the Alignment object #
	########################
	$gl = join("", @GL);
	$el = join("", @EL);
	$as = join("", @AS);
#print "$gl\n$as\n$el\n";
	my $ggaps = $gl =~ tr/-/-/;
	my $egaps = $el =~ tr/-/-/;
	my $match = $as;
	$match =~ s/[^|]//g;
	$match = length $match;
	my $positive = "";
	my $length = 0;
	my $alignment = EGlite::Alignment::new( $this->{SCORE}, $match, $positive, $length, 
		$gb, $ge, $eb, $ee, $gl, $el, $as, $ggaps, $egaps);
	my $alignmentbuf = EGlite::Alignment::new( $this->{SCORE}, $match, $positive, $length, 
		$gb, $ge, $eb, $ee, $gl, $el, $as, $ggaps, $egaps);
	$this->{ALIGNMENT} = $alignment;
	my $locus;
	$this->resetElements;
	my $elements = $this->{ELEMENTS};
	while (my $element = $this->nextElement) {
		if ($element->type eq "EXON") {
			$locus = index( $alignmentbuf->ea, ".");
			if ($locus == -1) { $locus = length $alignmentbuf->ea;}
		} elsif ($element->type eq "INTRON") {
			$alignmentbuf->ea =~ /\.+(\S+)/;
			my $next_exon = $1;
			if (defined $next_exon){
			    $locus = index( $alignmentbuf->ea, $next_exon);
			}
			if ($locus == -1) { warn "Possible error in nextElement()\n"; }
		} else {die "Unexpected element\n";}

#		print $alignmentbuf->ea."\n";
#		print $alignmentbuf->as."\n";
#		print $alignmentbuf->ga."\n";
		my $ea = substr( $alignmentbuf->ea, 0, $locus);
		my $as = substr( $alignmentbuf->as, 0, $locus);
		my $ga = substr( $alignmentbuf->ga, 0, $locus);
		$alignmentbuf->setEA( substr( $alignmentbuf->ea, $locus));
		$alignmentbuf->setAS( substr( $alignmentbuf->as, $locus));
		$alignmentbuf->setGA( substr( $alignmentbuf->ga, $locus));
		$element->setEA( $ea);
		$element->setAS( $as);
		$element->setGA( $ga);
#print "alignment:\n$ga\n$as\n$ea\n";
		shift @$elements;
		push @$elements, $element;
	}
	$this->resetElements;
	return 1;
}

sub _fastForward {
	my ($this) = @_;
	return 0 if $this->{REPORT_DONE};
	return 1 if $this->{LASTLINE} =~ /^>/;
	if ($this->{LASTLINE} =~ /^Parameters|^\s+Database:/) {
		$this->{REPORT_DONE} = 1;
		return 1;
	}
	my $FH = $this->{FH};
	while(<$FH>) {
		die if /^BLAST/;
		if ($_ =~ /^>|^Parameters|^\s+Database:/) {
			$this->{LASTLINE} = $_;
			$this->{REPORT_DONE} = 1;
			return 1;
		}
	}
	warn "Possible parse error in _fastForward in EGlite.pm\n";
}
sub _overload {
	my ($this) = @_;
	return $this->{GENOME} . " vs. " . $this->{EST};
}

###############################################################################
# EGlite::Element
###############################################################################
package EGlite::Element;
use overload '""' => 'name';
sub new {
	my $invocant = shift;
	my $class = ref($invocant) || $invocant;
	my $self = { @_ };          # Remaining args become attributes
	bless($self, $class);       # Bestow objecthood
	$self->{HSP_ALL_PARSED} = 0;
	return $self;
}
#	($self->{GENOME}, $self->{EST}, $self->{GB}, $self->{GE}, $self->{EB}, $self->{EE}, $self->{SCORE}, $self->{PI},
#		$self->{PARENT}) = @_;
sub name {
	my $self = shift;
	return $self->{GENOME}." vs ".$self->{PARENT}->{EST};
}
sub length {
	my $self = shift;
	return $self->{GE} - $self->{GB} + 1;
}
sub type            {shift->{TYPE}}
sub genome          {shift->{GENOME}}
sub est             {shift->{PARENT}->{EST}}
sub genomeBegin     {shift->{GB}}
sub genomeEnd       {shift->{GE}}
sub estBegin        {shift->{EB}}
sub estEnd          {shift->{EE}}
sub genomeAlignment {shift->{GL}}
sub estAlignment    {shift->{EL}}
sub alignmentString {shift->{AS}}
sub direction       {shift->{DIRECTION}}
sub gb              {shift->{GB}}
sub ge              {shift->{GE}}
sub eb              {shift->{EB}}
sub ee              {shift->{EE}}
sub ga              {shift->{GL}}
sub ea              {shift->{EL}}
sub as              {shift->{AS}}
sub score           {shift->{SCORE}}
sub pi              {shift->{PI}}
sub setEA           {
	my $self    = shift;
	$self->{EL} = shift;
}
sub setAS           {
	my $self    = shift;
	$self->{AS} = shift;
}
sub setGA           {
	my $self    = shift;
	$self->{GL} = shift;
}
sub alignment       {
	my $self    = shift;
	return $self->{GL}."\n".$self->{AS}."\n".$self->{EL};
}

###############################################################################
# EGlite::Exon
###############################################################################
package EGlite::Exon;
our @ISA = "EGlite::Element";
sub new {
	my $self = EGlite::Element->new(@_, TYPE => "EXON");
#($self->{GENOME}, $self->{EST}, $self->{GB}, $self->{GE}, $self->{EB}, $self->{EE}, $self->{SCORE}, $self->{PI},
#		$self->{PARENT}, "EXON") = @_;
	$self->{HSP_ALL_PARSED} = 0;
	return $self;
}

###############################################################################
# EGlite::Intron
###############################################################################
package EGlite::Intron;
our @ISA = "EGlite::Element";
sub new {
	my $self = EGlite::Element->new(@_, TYPE => "INTRON");
	$self->{HSP_ALL_PARSED} = 0;
	return $self;
}
sub est  {
	my $self = shift;
	return $self->{PARENT}->{EST};
}

###############################################################################
# EGlite::Alignment
###############################################################################
package EGlite::Alignment;
use overload '""' => '_overload';
sub new {
	my $alignment = bless {};
	($alignment->{SCORE}, 
		$alignment->{MATCH}, $alignment->{POSITIVE}, $alignment->{LENGTH},
		$alignment->{GB}, $alignment->{GE}, $alignment->{EB}, $alignment->{EE},
		$alignment->{GL}, $alignment->{EL}, $alignment->{AS}, $alignment->{GG}, $alignment->{EG}) = @_;
	$alignment->{PERCENT} = 0; #int(1000 * $alignment->{MATCH}/$alignment->{LENGTH})/10;
	return $alignment;
}
sub _overload {
	my $alignment = shift;
	return $alignment->ga."\n".$alignment->as."\n".$alignment->ea;
}
sub score           {shift->{SCORE}}
sub percent         {shift->{PERCENT}}
sub match           {shift->{MATCH}}
sub length          {shift->{LENGTH}}
sub genomeBegin     {shift->{GB}}
sub genomeEnd       {shift->{GE}}
sub estBegin        {shift->{EB}}
sub estEnd          {shift->{EE}}
sub genomeAlignment {shift->{GL}}
sub estAlignment    {shift->{EL}}
sub alignmentString {shift->{AS}}
sub genomeGaps      {shift->{GG}}
sub estGaps         {shift->{EG}}
sub gb              {shift->{GB}}
sub ge              {shift->{GE}}
sub eb              {shift->{EB}}
sub ee              {shift->{EE}}
sub ga              {shift->{GL}}
sub ea              {shift->{EL}}
sub as              {shift->{AS}}
sub gg              {shift->{GG}}
sub eg              {shift->{EG}}
sub setEA           {
	my $self    = shift;
	$self->{EL} = shift;
}
sub setAS           {
	my $self    = shift;
	$self->{AS} = shift;
}
sub setGA           {
	my $self    = shift;
	$self->{GL} = shift;
}

###############################################################################
# EGlite::Multi
###############################################################################
package EGlite::Multi;
sub new {
	my ($class, $fh) = @_;
	if (ref $fh !~ /GLOB/)
		{die "EGlite error: new expects a GLOB reference not $fh\n"}
	my $this = bless {};
	$this->{FH} = $fh;
	return $this;
}
sub nextReport {
	my ($this) = @_;
	my $FH = $this->{FH};
	my $blast = new EGlite($this->{FH});
	return $blast;
}

1;
__END__

=head1 NAME

EGlite - Lightweight est2genome parser

=head1 SYNOPSIS

 use EGlite;
 
 # single E2G report
 
 my $report = new EGlite(\*STDIN);
 $report->genome;
 $report->est;
 $report->score;
 $report->pi;
 $report->alignment;
 $report->genomeBegin;
 $report->genomeEnd;
 $report->estBegin;
 $report->estEnd;
 $report->genomeDirection;
 $report->estDirection;
 $report->gb;
 $report->ge;
 $report->eb;
 $report->ee;
 $report->hasIntron;
 
 if( $report->hasIntron) {
     #do something with introns
     if (my $alignment = $report->alignment) {
         $alignment->score;
         $alignment->percent;
         $alignment->match;
         $alignment->length;
         $alignment->genomeBegin;
         $alignment->genomeEnd;
         $alignment->estBegin;
         $alignment->estEnd;
         $alignment->genomeAlignment;
         $alignment->estAlignment;
         $alignment->alignmentString;
         $alignment->genomeGaps;
         $alignment->estGaps;
     }
 }
 
 while(my $element = $report->nextElement) {
     $element->name;
     $element->length;
     $element->type; #"EXON" or "INTRON"
     $element->genome;
     $element->est;
     $element->genomeBegin;
     $element->genomeEnd;
     $element->estBegin;
     $element->estEnd;
     $element->genomeAlignment;
     $element->estAlignment;
     $element->alignmentString;
     $element->gb;
     $element->ge;
     $element->eb;
     $element->ee;
     $element->ga;
     $element->ea;
     $element->as;
     $element->score;
     $element->pi;
 }

 $report->resetElements; 
 while($report->hasMoreElements) {
     my $element = $report->nextElement;
 }

 $report->resetElements;
 my $count = $report->countElements;
 for( my $i=0; $i < $count; $i++){
     my $element = $report->nextElement;
 }
 
 # multiple (concatenated) E2G reports
 
 my $multiple_report = new EGlite::Multi(\*STDIN);
 while(my $report = $multiple_report->nextReport) {
     while(my $element = $report->nextElement) {
         #whatever
     }
 }


=head1 DESCRIPTION

EGlite is a (pseudo-)lightweight package for parsing est2genome reports. 
est2genome is a program to align Expressed Sequence Tags (EST) against a 
genomic sequence. The output gives detailed positional information about
(potential) exons and introns supported by the EST sequence. EGlite parses
this output and provides a handle on the relevant and required information. 

=head2 Object

EGlite has three kinds of objects, the report, the element, and the alignment. 
The element has two subtypes: exon and intron. To create a new report, you 
pass a filehandle reference to the EGlite constructor.

 my $report = new EGlite(\*STDIN); # or any other filehandle

The report has many attributes and methods.

 $report->genome;           # genomic sequence name 
 $report->est;              # EST sequence name
 $report->score;            # Score of the overall alignment
 $report->pi;               # Percentage identity in the overall alignment
 $report->alignment;        # The alignment object
 $report->genomeBegin;      # The start co-ordinate of the genomic alignment
 $report->genomeEnd;        # The end co-ordinate of the genomic alignment
 $report->estBegin;         # The start co-ordinate of the EST alignment
 $report->estEnd;           # The end co-ordinate of the EST alignment
 $report->genomeDirection;
 $report->estDirection;
 $report->gb;
 $report->ge;
 $report->eb;
 $report->ee;
 $report->hasIntron;        # Does the alignment contain an intron?
 $report->resetElements;    # Resets the pointer of current element to the first element
 $report->hasMoreElements;  # Does the alignment has more elements?
 $report->nextElement;      # Next element in the alignment

 my $alignment = $report->alignment;
 print "$alignment"; 

 The output of such code might look like this:

 agcag-t
 | ||  | 
 atcatct

 The traditional access through a lightweight parser is using:

 while(my $element = $report->nextElement) {
     #do something with $element
 }

 Alternatively,

 $report->resetElements; 
 while($report->hasMoreElements) {
     my $element = $report->nextElement;
 }

 or

 $report->resetElements;
 my $count = $report->countElements;
 for( my $i=0; $i < $count; $i++){
     my $element = $report->nextElement;
 }

 give you more control over the iterations. You can stop the iteration
 at your convenience using these two constructs.

=head2 Concatenated est2genome reports

You can step through multiple est2genome reports if you have a file of concatenated
est2genome reports using the following construct.

 my $multiple_report = new EGlite::Multi(\*STDIN);
 while(my $report = $multiple_report->nextReport) {
     while(my $element = $report->nextElement) {
     }
 }


=head1 AUTHOR

Manimozhiyan Arumugam (mozhiyan@cse.wustl.edu)

=head1 ACKNOWLEDGEMENTS

This program was adapted from Ian Korf's BPlite.pm. 
This software was developed at the Computational Genomics Laboratory at Washington
Univeristy, St. Louis, MO.

=head1 DISCLAIMER

This software is provided "as is" without warranty of any kind.

=cut










