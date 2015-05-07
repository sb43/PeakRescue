package PeakRescue::GetPeak; 
use PeakRescue;
our $VERSION = PeakRescue->VERSION;

use strict;
use Bio::DB::Sam;
use Bio::DB::Sam::Constants;



use FindBin qw($Bin);
use List::Util qw(max min);
use File::Path qw(mkpath remove_tree);
use Capture::Tiny qw(:all);
use File::Basename;
use File::Spec;
use Data::Dumper;
use Const::Fast qw(const);
use Log::Log4perl;
Log::Log4perl->init("$Bin/../config/log4perl.gt.conf");
my $log = Log::Log4perl->get_logger(__PACKAGE__);

use PeakRescue::Base;

const my $PROPER_PAIRED => 0x2;
const my $UNMAPPED => 0x4;
const my $REVERSE_STRAND => 0x10;
const my $MATE_REVERSE_STRAND => 0x20;
const my $NOT_PRIMARY_ALIGN => 0x100;
const my $MATE_UNMAPPED => 0x0008;
const my $READ_PAIRED => 0x0001;
const my $FIRST_IN_PAIR => 0x40;
const my $MAX_PILEUP_DEPTH => '1000000';
const my $SUPP_ALIGNMENT => 0x800;
const my $DUP_READ => 0x400;
const my $VENDER_FAIL => 0x200;



sub new {
	my ($class,$options)=@_;
	$log->trace('new');
	my $self={};
	bless $self, $class;
	$self->_init($options);
	$self->_get_gene_interval();
	#$self->_get_max_coverage();
	return $self;
}

=head2 _init
populate object with  necessary file names
Inputs
=over 2
=item options -user input paramaters
=back
=cut

sub _init {
	my ($self,$options)=@_;
	my ($file_name,$dir_name,$suffix) = fileparse($options->{'bam'},qr/\.[^.]*/);
	$options->{'s'}=$suffix;
	$options->{'d'}=$dir_name;
  if ($options->{'o'} && ! (-d $options->{'o'})) {
  	$log->debug("Creating dir:".$options->{'o'});
		mkpath($options->{'o'});
  }
	elsif(!$options->{'o'}) {
   $options->{'o'}=$dir_name;
  }
  elsif (!(-d $options->{'o'})) {
  	$log->logcroak("Unable to create directory");
  }
	$options->{'f'}=$options->{'o'}.$file_name;	
	mkpath($options->{'o'}.'/'.'tmp_peak');
	$options->{'tmpdir'}=$options->{'o'}.'/'.'tmp_peak';
	$self->{'options'} = $options;
	
	$log->debug("Using bam file:".$options->{'bam'});
	$log->debug("Using bed file:".$options->{'bed'});
	return;
}

sub _get_bam_object {
	my ($self,$bam)=@_;
		my $bam = Bio::DB::Sam->new(-bam => $bam,
															-fasta => $self->options->{'g'},
															-expand_flags => 1,
															-split_splices => 1
															);
		$bam->max_pileup_cnt($MAX_PILEUP_DEPTH); # chanages default coverage cap from 8000 to $new_value
		return $bam;
}



sub _get_gene_interval { 
	my ($self)=@_;
	my $peak_data;
	my($bam_object)=$self->_get_bam_object($self->options->{'bam'});
	my $bam_header=$bam_object->header;
	my $counter;
	open(my $bed_fh, '<' , $self->options->{'gb'}) || $log->logcroak ("unable to open bed file $!0");
	while (<$bed_fh>) {
		chomp;
		my ($chr,$start,$end,$gene)=(split "\t", $_) [0,1,2,3];
		$counter++;
		if($counter % 1000 == 0) {
			$log->debug("Calculated peak for $gene:==>".$counter);
		}
		my($gene_bam_object) = $self->_get_gene_bam_object($bam_object,$chr,$start,$end,$gene,$bam_header);
		my ($coverage) = $gene_bam_object->features(-type => 'coverage', -seq_id => $chr, -start => $start, -end => $end);
		my @data  = $coverage->coverage;
		$peak_data->{$gene}=max(@data);
		undef $gene_bam_object;
	}
	$log->debug("Completed peak calculation===>".$counter.' gene');
	open(my $tmp_fh, ">tmp_peak.txt");
	$self->print_peak($peak_data,$tmp_fh);
	close($tmp_fh);
	
}


=head2 _get_max_coverage
get maximum coverage base in an interval

Inputs
=over 2
=item sam_object - Bio::DB sam object
=back
=cut

sub _get_max_coverage {
	my($self)=@_;
	my ($peak_data,$cov_array,$chr,$start,$end,$gene,$counter,$tmp_gene);
	my($bam_object)=$self->_get_bam_object();
	if(	$bam_object) {$log->debug("Bam object created");}
	my $bam_header=$bam_object->header->text;
	open(my $bed_fh, '<' , $self->options->{'bed'}) || $log->logcroak ("unable to open bed file $!0");	
	open(my $tmp_fh, ">tmp_peak.txt");
	while (<$bed_fh>) {
		chomp;
		($chr,$start,$end,$gene)=(split "\t", $_) [0,1,2,3];
		my($gene_bam) = $self->_get_gene_bam($bam_object,$chr,$start,$end,$gene,$bam_header,$tmp_fh);
		# Currently filter is ignored. In reality, we should
    # turn filter into a callback and invoke it on each 
    # position in the pileup.
    
		my ($coverage) = $bam_object->features(-type => 'coverage', -seq_id => $chr, -start => $start, -end => $end);
		my @data  = $coverage->coverage;
		if(defined $peak_data->{$gene}){
			push(@$cov_array,max(@data));
		}
		else{
			if(defined $peak_data->{$tmp_gene}) {
				$counter++;
				$peak_data->{$tmp_gene}=$cov_array;
				#print hash after every 1000 genes
				if($counter % 1000 == 0) {
					$self->print_peak($peak_data,$tmp_fh);
					$peak_data=();
					$log->debug("Completed peak calculation for: $counter genes");
				}
			}
			$cov_array=();
			push(@$cov_array,max(@data));
			$peak_data->{$gene}= $cov_array;
			$tmp_gene=$gene;
		}
	}
			$peak_data->{$gene}= $cov_array;
			$self->print_peak($peak_data,$tmp_fh);
			$log->debug("Completed peak calculation for total:".($counter+1)."genes");
			close($tmp_fh);
			$peak_data=();
return;
}

sub print_peak {
	my ($self,$hash,$fh)=@_;
	foreach my $key (sort keys %$hash) {
		if(defined($hash->{$key}))
		{	
			print $fh "$key\t".$hash->{$key}."\n";	
		}
	}
}



sub _get_gene_bam_object {
	my ($self,$bam_object,$chr,$start,$end,$gene,$header)=@_;
	my $tmp_gene_file=$self->options->{'tmpdir'}.'/tmp_gene.bam';
	#my ($fh)=PeakRescue::Base::_create_fh([$tmp_gene_file],1);
	#my $fh_str=@$fh[0];
	#close($fh_str);
	my $mode = 'w';
	my $bam = Bio::DB::Bam->open($tmp_gene_file,$mode);
	$bam->header_write($header);
	my @alignments = $bam_object->get_features_by_location(-seq_id => $chr,
																									 -start  => $start,
																									 -end    => $end);
																									 																							 
	 for my $align (@alignments) {
			my $str = $align->aux;
			if ($str =~m/XF:Z:$gene/) {
			  	# alignments filtered  by given gene 
					# these objects are B:D:AlignWrapper objects, not B:D:Alignment 
					# objects, therefore I can't simply just write them out to the 
					# the low level bam object
					# there is no official way to convert from AlignWrapper to 
					# Alignment (even though AlignWrapper wraps around Alignment).
					
					# therefore, I have to manually extract the Alignment object.
					# Using Dumper, identified the Alignment object as the value 
					# under the object's key 'align'
					# and then write that to the output bam file [ as mentioned in biotoolbox - split_bam_by_isize.pl script ]
			  $bam->write1($align->{'align'});
			}
	 }
	 # this is required to put EOF marker before creating bam object
	undef $bam;
	
	#system("samtools index $tmp_gene_file");
	Bio::DB::Bam->index_build($tmp_gene_file);
	my ($gene_bam)=$self->_get_bam_object($tmp_gene_file);
	return $gene_bam;
}










sub options {
shift->{'options'};
}




1;