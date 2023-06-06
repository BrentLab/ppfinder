#! perl -w
package Sim4Parser;
use strict;
use overload '""' => '_overload';
###############################################################################
# Sim4Parser
###############################################################################
sub new {
	my ($class, $fh) = @_;
	if (ref $fh !~ /GLOB/)
		{die "Sim4Parser error: new expects a GLOB reference not $fh\n"}
	my $this = bless {};
	$this->{FH} = $fh;
	$this->{LASTLINE} = "";

	# Parse the header
	if ($this->_parseHeader) {
		# there are alignments
		$this->{REPORT_DONE} = 0;
		if(not defined $this->{GENOME}) {die "Incorrect Sim4 output format: Choose A=4\n";}
		if(not defined $this->{EST}) {die "Incorrect Sim4 output format: Choose A=4\n";}

		# Parse the alignment
		my $FH = $this->{FH};
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
#sub genomeDirection {return "N/A"}
#sub estDirection    {return "N/A"}
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
	my $counter = 0;
	my @genome_array = ();
	my @est_array    = ();
	my %intron_orientation = ('->' => 0, '<-' => 0, '==' => 0);
	$this->{GEN_DIRECTION} = '+';
	my ($last_exon_gb, $last_exon_ge, $last_exon_eb, $last_exon_ee, $last_exon_dir);

	while(<$FH>) {

		# File name information
		if ($_ =~ /seq1 = (.*), ([0-9]+) bp/ ) {
			# Can do something if we need the filenames at all!	
		}
		elsif ($_ =~ /seq2 = (.*), ([0-9]+) bp/ ) {
		}

		# Is the genomic sequence reverse complemented?
		elsif ($_ =~ /\(complement\)/) {
			$this->{GEN_DIRECTION} = '-';
		}

		# Fasta header information
		elsif ($_ =~ /^>(.*)$/ ){
			if($counter == 0){
				$this->{EST} = $1;
				$counter++;
			}elsif($counter == 1){
				$this->{GENOME} = $1;
			}else{
				die "3rd occurrence of > at start of line. Quitting...\n";
			}
		}

		# Alignment coordinates information
		elsif ($_ =~ /([0-9]+)-([0-9]+)  \(([0-9]+)-([0-9]+)\)   ([0-9]+)% (->|<-|==|--)/ ) {
			my ($pi, $gb, $ge, $eb, $ee, $direction) = ($5, $3, $4, $1, $2, $6);	
			($intron_orientation{$direction})++;
			if (defined $last_exon_gb) {
				my $intron = Sim4Parser::Intron::new( GENOME=>$this->{GENOME}, 
								EST=>$this->{EST}, 
								GB=>($last_exon_ge + 1), 
								'GE'=>($gb - 1), 
								EB=>($last_exon_ee + 1), 
								'EE'=>($eb - 1), 
								SCORE=>"N/A", 
								PI=>0, 
								DIRECTION=>$last_exon_dir,
								PARENT=>$this);
				push @introns, $intron;
				push @elements, $intron;
			}
			my $exon = Sim4Parser::Exon::new( GENOME=>$this->{GENOME}, 
							EST=>$this->{EST}, 
							GB=>$gb, 
							'GE'=>$ge, 
							EB=>$eb, 
							EE=>$ee, 
							SCORE=>"N/A", 
							PI=>$pi, 
							DIRECTION=>$direction,
							PARENT=>$this);
			push @exons, $exon;
			push @elements, $exon;
			push @est_array,    ($eb, $ee);
			push @genome_array, ($gb, $ge);
			($last_exon_gb, $last_exon_ge, $last_exon_eb, $last_exon_ee, $last_exon_dir) = ($gb, $ge, $eb, $ee, $direction);
		}

		# Last exon information
		elsif ($_ =~ /([0-9]+)-([0-9]+)  \(([0-9]+)-([0-9]+)\)   ([0-9]+)%/ ) {
			my ($pi, $gb, $ge, $eb, $ee) = ($5, $3, $4, $1, $2);	
			if (defined $last_exon_gb) {
				my $intron = Sim4Parser::Intron::new( GENOME=>$this->{GENOME}, 
								EST=>$this->{EST}, 
								GB=>($last_exon_ge + 1), 
								'GE'=>($gb - 1), 
								EB=>($last_exon_ee + 1), 
								'EE'=>($eb - 1), 
								SCORE=>"N/A", 
								PI=>0, 
								DIRECTION=>$last_exon_dir,
								PARENT=>$this);
				push @introns, $intron;
				push @elements, $intron;
			}
			my $exon = Sim4Parser::Exon::new( GENOME=>$this->{GENOME}, 
							EST=>$this->{EST}, 
							GB=>$gb, 
							'GE'=>$ge, 
							EB=>$eb, 
							EE=>$ee, 
							SCORE=>"N/A", 
							PI=>$pi, 
							DIRECTION=>"N/A",
							PARENT=>$this);
			push @exons, $exon;
			push @elements, $exon;
			push @est_array,    ($eb, $ee);
			push @genome_array, ($gb, $ge);
			$this->{LASTLINE} = $_; 
			$this->{SCORE}    = "N/A";
			$this->{EXONS}    = \@exons;
			$this->{INTRONS}  = \@introns;
			$this->{ELEMENTS} = \@elements;
			$this->resetElements; # Sets ELEMENTBUFFER
			($this->{EST_BEGIN},    $this->{EST_END})    = ($est_array[0],    $est_array[$#est_array]);
			($this->{GENOME_BEGIN}, $this->{GENOME_END}) = ($genome_array[0], $genome_array[$#genome_array]);
			if ($intron_orientation{'->'} > $intron_orientation{'<-'}) {
				$this->{SPLICE_DIRECTION} = '+';
			} elsif ($intron_orientation{'->'} < $intron_orientation{'<-'}) {
				$this->{SPLICE_DIRECTION} = '-';
			} else {
				$this->{SPLICE_DIRECTION} = '.';
			}
			$this->{EST_DIRECTION} = '+';
			return 1;
		}

		# Fix this soon
		elsif ($_ =~ /Note Best alignment is between (reversed|forward) est and (reversed|forward) genome, (but|and) splice\s+sites imply\s+(forward gene|REVERSED GENE)/ ) {
		}
	}
	return 0 unless (defined $this->{GENOME} && defined $this->{EST});
}

sub _parseAlignment {
	my ($this) = @_;
	return 0 if $this->{ALIGN_ALL_PARSED};

	my $est = $this->{EST};
	my $genome = $this->{GENOME};
	my $flag = 1;
	my $FH  = $this->{FH};
	
	#######################
	# get alignment lines #
	#######################
	my @algnline;
	#$_ = $this->{LASTLINE};
	while (<$FH>) {
		if ($_ =~ /^\s*(0)[\s\.:]+$/) {
			my $l0 = <$FH>; push @algnline, $l0;  # grab/store the genome line
			my $l1 = <$FH>; push @algnline, $l1; # grab/store the alignment line
			my $l2 = <$FH>; push @algnline, $l2; # grab/store the est line
#print "$l0$l1$l2";
			last;
		}
	}
	while(<$FH>) {
		die if /^BLAST/;
#print $_;

		if ($_ !~ /\S/)            {next}

		# Coordinate scale line starts with 0 for a new HSP
		elsif ($_ =~ /^\s*(0)[\s\.:]+$/) {
			push @algnline, "=========";
			push @algnline, "=========";
			push @algnline, "=========";
			my $l0 = <$FH>; push @algnline, $l0;  # grab/store the genome line
			my $l1 = <$FH>; push @algnline, $l1; # grab/store the alignment line
			my $l2 = <$FH>; push @algnline, $l2; # grab/store the est line
#print "$l0$l1$l2";
		}

		# Coordinate scale line
		elsif ($_ =~ /^\s*([1-9]*[05]*0)[\s\.:]+$/) {
			my $l0 = <$FH>; push @algnline, $l0;  # grab/store the genome line
			my $l1 = <$FH>; push @algnline, $l1; # grab/store the alignment line
			my $l2 = <$FH>; push @algnline, $l2; # grab/store the est line
#print "$l0$l1$l2";
		}

		# End of this report!
		elsif ($_ =~ /NOTE: sim4 report stops here/ ) {
			$this->{LASTLINE} = $_;
#print "END\n";
			last;
		}
		$this->{LASTLINE} = $_;
	}
	$this->{ALIGN_ALL_PARSED} = 1;
	
	#########################
	# parse alignment lines #
	#########################
	my ($gl, $el, $as) = ("", "", "");
	my ($gb, $ge, $eb, $ee) = (0,0,0,0);
	my (@GL, @EL, @AS); # for better memory management
			
	for(my $i=0;$i<@algnline;$i+=3) {
		#warn $algnline[$i], $algnline[$i+2];

		# Parse the lines

		# Genomic line
		if ($algnline[$i] eq "=========") {
			push @GL, "========="; push @EL, "========="; push @AS, "=========";
			next;
		}
		$algnline[$i]   =~ /^\s+(\d+) (.+)$/;
		return 0 unless defined $2;
		$el = $2; $eb = $1 unless $eb;

		# Set the offset based on this line
		my $offset = index($algnline[$i], $el);

		$as = substr($algnline[$i+1], $offset, CORE::length($el))
			if $algnline[$i+1];
		
		$algnline[$i+2] =~ /^\s+(\d+) (.+)$/;
		return 0 unless defined $2;
		$gl = $2; $gb = $1 unless $gb;
		
		my $gap = -1;
		while (($gap = index($as, "-", $gap)) != -1) {
			if (substr($gl, $gap, 1) eq " ") {
				substr($gl, $gap, 1) = "-";
			} elsif (substr($el, $gap, 1) eq " ") {
				substr($el, $gap, 1) = "-";
			}
			substr($as, $gap, 1) = " ";
		}
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
	my $ggaps = $gl;
	my $egaps = $el;
	$ggaps = $ggaps =~ s/\s+/ /g;
	$egaps = $egaps =~ s/\s+/ /g;
	my $match = $as;
	$match =~ s/[^|]//g;
	$match = length $match;
	my $positive = "";
	my $length = abs($this->{GENOME_END} - $this->{GENOME_BEGIN}) + 1;

	#Calculate global PI
	my $as_stripped = $as;
	$as_stripped =~ s/[>\.<]+//g;
	my $as_length   = length $as_stripped;
	my $global_pi = int(1000*($match)/$as_length)/10;
	$this->{PI} = $global_pi;
	
	my $alignment = Sim4Parser::Alignment::new( $this->{SCORE}, $global_pi, $match, $positive, $length, 
		$this->{GENOME_BEGIN}, $this->{GENOME_END}, $this->{EST_BEGIN}, $this->{EST_END}, $gl, $el, $as, $ggaps, $egaps);
	my $alignmentbuf = Sim4Parser::Alignment::new( $this->{SCORE}, $global_pi, $match, $positive, $length, 
		$gb, $ge, $eb, $ee, $gl, $el, $as, $ggaps, $egaps);
	$this->{ALIGNMENT} = $alignment;

	my $locus;
	$this->resetElements;

	# Set individual element objects
	my $elements = $this->{ELEMENTS};
	while (my $element = $this->nextElement) {
		if ($element->type eq "EXON") {
			my ($min_locus, $intron_symbol) = (-1, "");
			foreach my $symbol (">", "<", "?", "=") {
				$locus = index($alignmentbuf->as, $symbol);
				if ($locus != -1 && ($locus < $min_locus || $min_locus == -1)) {
					$min_locus = $locus;
					$intron_symbol = $symbol;
				}
			}
			if ($min_locus == -1) { $min_locus = (length $alignmentbuf->ea);}

			$locus = $min_locus;
		} elsif ($element->type eq "INTRON") {
			$alignmentbuf->ga =~ /(\S{3}\.\.\.\S{3}\S+|\={9}\S+)/;
#print "BUF:".$alignmentbuf->ga."\n";

			my $next_exon = $1;
			if (!defined($next_exon)) {
				die "Cannot get the next exon";
			}
			# You do it this way since there might be repeats that will kill you. GGA...CAAAAATTTA will return the AAA after C instead of AAA after CAA. 
			$locus = index( $alignmentbuf->ga, $next_exon) + 9;
			if ($locus == -1) { warn "Possible error in nextElement()\n"; }
		} else {die "Unexpected element\n";}

#		print $alignmentbuf->ea."\n";
#		print $alignmentbuf->as."\n";
#		print $alignmentbuf->ga."\n";
		my $ea = substr( $alignmentbuf->ea, 0, $locus);
		my $as = substr( $alignmentbuf->as, 0, $locus);
		my $ga = substr( $alignmentbuf->ga, 0, $locus);
#print STDERR $alignmentbuf->ea."\t$locus\n";
		if( $locus == (length $alignmentbuf->ea) ) {
			$alignmentbuf->setEA( substr( $alignmentbuf->ea, $locus));
			$alignmentbuf->setAS( substr( $alignmentbuf->as, $locus));
			$alignmentbuf->setGA( substr( $alignmentbuf->ga, $locus));
		}else{
			$alignmentbuf->setEA( substr( $alignmentbuf->ea, $locus + 0));
			$alignmentbuf->setAS( substr( $alignmentbuf->as, $locus + 0));
			$alignmentbuf->setGA( substr( $alignmentbuf->ga, $locus + 0));
		}
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
	warn "Possible parse error in _fastForward in Sim4Parser.pm\n";
}

sub _overload {
	my ($this) = @_;
	return "\'".$this->{GENOME}."\' vs. \'".$this->{EST}."\'";
}

###############################################################################
# Sim4Parser::Element
###############################################################################
package Sim4Parser::Element;
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
	return "\'".$self->{GENOME}."\' vs \'".$self->{PARENT}->{EST}."\'";
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
sub strval {
	my $self = shift;
	my $string = sprintf "%9d %s %9d\n%9s %s %9s", $self->gb, $self->ga, $self->ge, " ", $self->as, " ";
	if ($self->type eq "EXON") {
		$string = sprintf "%s\n%9d %s %9d", $string, $self->eb, $self->ea, $self->ee;
	}
	return $string;
}

###############################################################################
# Sim4Parser::Exon
###############################################################################
package Sim4Parser::Exon;
our @ISA = "Sim4Parser::Element";
sub new {
	my $self = Sim4Parser::Element->new(@_, TYPE => "EXON");
#($self->{GENOME}, $self->{EST}, $self->{GB}, $self->{GE}, $self->{EB}, $self->{EE}, $self->{SCORE}, $self->{PI},
#		$self->{PARENT}, "EXON") = @_;
	$self->{HSP_ALL_PARSED} = 0;
	return $self;
}

###############################################################################
# Sim4Parser::Intron
###############################################################################
package Sim4Parser::Intron;
our @ISA = "Sim4Parser::Element";
sub new {
	my $self = Sim4Parser::Element->new(@_, TYPE => "INTRON");
	$self->{HSP_ALL_PARSED} = 0;
	return $self;
}
sub est  {
	my $self = shift;
	return $self->{PARENT}->{EST};
}

###############################################################################
# Sim4Parser::Alignment
###############################################################################
package Sim4Parser::Alignment;
use overload '""' => '_overload';
sub new {
	my $alignment = bless {};
	($alignment->{SCORE}, $alignment->{PI},
		$alignment->{MATCH}, $alignment->{POSITIVE}, $alignment->{LENGTH},
		$alignment->{GB}, $alignment->{GE}, $alignment->{EB}, $alignment->{EE},
		$alignment->{GL}, $alignment->{EL}, $alignment->{AS}, $alignment->{GG}, $alignment->{EG}) = @_;
	return $alignment;
}
sub _overload {
	my $alignment = shift;
	return $alignment->ga."\n".$alignment->as."\n".$alignment->ea;
}
sub score           {shift->{SCORE}}
sub pi              {shift->{PI}}
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
# Sim4Parser::Multi
###############################################################################
package Sim4Parser::Multi;
sub new {
	my ($class, $fh) = @_;
	if (ref $fh !~ /GLOB/)
		{die "Sim4Parser error: new expects a GLOB reference not $fh\n"}
	my $this = bless {};
	$this->{FH} = $fh;
	return $this;
}
sub nextReport {
	my ($this) = @_;
	my $FH = $this->{FH};
	my $blast = new Sim4Parser($this->{FH});
	return $blast;
}

1;
__END__

=head1 NAME

Sim4Parser - Lightweight est2genome parser

=head1 SYNOPSIS

 use Sim4Parser;
 
 # single E2G report
 
 my $report = new Sim4Parser(\*STDIN);
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
         $alignment->pi;
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
 
 my $multiple_report = new Sim4Parser::Multi(\*STDIN);
 while(my $report = $multiple_report->nextReport) {
     while(my $element = $report->nextElement) {
         #whatever
     }
 }


=head1 DESCRIPTION

Sim4Parser is a (pseudo-)lightweight package for parsing est2genome reports. 
est2genome is a program to align Expressed Sequence Tags (EST) against a 
genomic sequence. The output gives detailed positional information about
(potential) exons and introns supported by the EST sequence. Sim4Parser parses
this output and provides a handle on the relevant and required information. 

=head2 Object

Sim4Parser has three kinds of objects, the report, the element, and the alignment. 
The element has two subtypes: exon and intron. To create a new report, you 
pass a filehandle reference to the Sim4Parser constructor.

 my $report = new Sim4Parser(\*STDIN); # or any other filehandle

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

 my $multiple_report = new Sim4Parser::Multi(\*STDIN);
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

