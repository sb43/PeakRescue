package PeakRescue::GetPeak; 
use PeakRescue;
our $VERSION = PeakRescue->VERSION;

use strict;
use Bio::DB::Sam;

use Tabix;
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
	$self->_do_max_peak_calculation();
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
	$options->{'f'}=$file_name;
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
	mkpath($options->{'o'}.'/'.'tmp_peak');
	$options->{'tmpdir'}=$options->{'o'}.'/'.'tmp_peak';
	$self->{'options'} = $options;
	
	$log->debug("Using bam file:".$options->{'bam'});
	$log->debug("Using bed file:".$options->{'bed'});
	$log->debug("Using temp dir:".$options->{'tmpdir'});
	return;
}


=head2  _get_bam_object
create biodbsam object for a given bam file
Inputs
=item bam_file - bam file to create object
=over 2
=back
=cut

sub _get_bam_object {
	my ($self,$bam_file)=@_;
		my $bam = Bio::DB::Sam->new(-bam => $bam_file,
															-fasta => $self->options->{'g'},
															-expand_flags => 1,
															-split_splices => 1
															);
		$bam->max_pileup_cnt($MAX_PILEUP_DEPTH); # changes default coverage cap from 8000 to $new_value
		return $bam;
}

=head2  _do_max_peak_calculation
main sub which does peak calculation using other subs
Also stores progress for every peak
Inputs
=over 2
=back
=cut

sub _do_max_peak_calculation { 
	my ($self)=@_;
	my($bam_object)=$self->_get_bam_object($self->options->{'bam'});
	my($tabix)=$self->_get_tabix_object();
	open(my $peak_fh, '>>', $self->options->{'o'}.'/'.$self->options->{'f'}.'_peak.txt');
	# check progress required as this is long process if beaks in between can be restarted from where it left using progress file
	my $progress_file=$self->options->{'o'}.'/'.'progress.txt';
	if( -e $progress_file) { $log->info("Progress file exists, analysis will be resumed, to restart analysis please remove progress file"); }
	open(my $progress_fh, '>>', $progress_file);
	open(my $read_progress, '<', $progress_file);
	my @chr_analysed=<$read_progress>; 
	close($read_progress);
	
	# end of progress check	
	$log->logcroak("Unable to create tabix object") if (!$tabix);
	# get chromosome names
	my @chrnames=$tabix->getnames;
	###
	foreach my $chr (@chrnames) {
	next if (grep (/^$chr$/, @chr_analysed));
	$log->debug(">>>>>>>Calculating peak for chromosome: $chr ");
	my ($peak_data, $counter);
	my $res = $tabix->query($chr);
		while(my $record = $tabix->read($res)){		
			chomp;
			my ($chr,$start,$end,$gene)=(split "\t", $record) [0,1,2,3];
			$counter++;
			if($counter % 1000 == 0) {
				$log->debug("Calculated peak for chromosome: $chr ==>".$counter.' genes');
				# can be used if hash is running out of memory
				#$self->_print_peak($peak_data,$peak_fh);
				#$peak_data=();
			}
			my($max_peak)=$self->_get_gene_bam_object($bam_object,$chr,$start,$end,$gene);
			$peak_data->{$gene}=$max_peak;
		}
	# print peak data for each chromosome
	$self->_print_peak($peak_data,$peak_fh);
	$log->debug(">>>>>>>>> Completed peak calculation for chromosome: $chr ===>".$counter.' gene >>>>>>>>>');
	print $progress_fh $chr."\n";
}
	close($peak_fh);
	close($progress_fh);
	
	#PeakRescue::Base::cleanup_dir($self->options->{'tmpdir'});
	
}

=head2  _get_tabix_object
create tabix object for a user provided bed file
Inputs
=over 2
=item bed
=back
=cut

sub _get_tabix_object {
	my ($self)=shift;
	my $tabix_obj;
	my $bed_file = $self->options->{'bed'};
	if (! -e $bed_file) {
		$log->logcroak("Unable to find file : $bed_file ");
	}
	my $tmp_bed = $self->options->{'tmpdir'}.'/tmpbed_sorted.bed';
	my ($out,$stderr,$exit) = capture{system("bedtools sort -i $bed_file | bgzip >$tmp_bed.gz && tabix -p bed $tmp_bed.gz ")};	
		if ($exit) {
			$log->logcroak("Unable to create tabix bed $stderr");
		}
		else {
			$tabix_obj = new Tabix(-data => "$tmp_bed.gz");
			$log->debug("Tabix object created successfully");
			return $tabix_obj;
		}
	return $tabix_obj;
}

=head2  _get_gene_bam_object
fetch reads
Inputs
=over 2
=item bam_object - Bio::DB bam object
=item region - chr:start-end format region to get reads
=item gene - gene name 
=back
=cut

sub _get_gene_bam_object {
my ($self,$bam_object,$chr,$start,$end,$gene)=@_;
	my $tmp_gene_file=$self->options->{'tmpdir'}.'/tmp_gene.bam';
	# create bio db bam object
	my $bam = Bio::DB::Bam->open($tmp_gene_file,'w');
	
	#write header
	$bam->header_write($bam_object->header);
	# create temp file
	my $read_flag=undef;
	$bam_object->fetch("$chr:$start-$end", sub {
		my $a = shift;
		my $flags = $a->flag;
		return if $flags & $NOT_PRIMARY_ALIGN;
		return if $flags & $VENDER_FAIL;
		return if $flags & $UNMAPPED;
		return if $flags & $DUP_READ;
		#return if $flags & $SUPP_ALIGNMENT;
		# exact match with XF tag value
		if ($a->get_tag_values('XF') eq $gene) {
			$bam->write1($a->{'align'});
			$read_flag=1;
		}elsif($a->aux=~ m/XF:Z:$gene/){
			$bam->write1($a->{'align'});
			$read_flag=1;
		}
	});
	undef $bam;
	# read flag is undefined 
	if (!$read_flag) {
		#print "No reads for :  $gene\n";
		return "0";
	} 
	Bio::DB::Bam->index_build($tmp_gene_file);
	# create tmp file
	my $tmp_gene_sorted=$self->options->{'tmpdir'}.'/tmp_gene_name_sorted';
	# name sorted bam required for clipOver
	Bio::DB::Bam->sort_core(1,$tmp_gene_file,$tmp_gene_sorted);
	my($clipped_bam)=$self->_run_clipOver("$tmp_gene_sorted.bam");
	my ($gene_bam)=$self->_get_bam_object($clipped_bam);	
	my ($coverage) = $gene_bam->features(-type => 'coverage', -seq_id => $chr, -start => $start, -end => $end);
	return max($coverage->coverage);
	
	
	#option2 samtools mpipeup and GATK [run separately] too slow.
  # need samtools 1.1 or above which takes care of overlapping read pairs... 
	#my $cmd= "~/software/samtools-1.2/samtools mpileup $tmp_gene_file -d $MAX_PILEUP_DEPTH -A -f ".$self->options->{'g'}. " -r $region --no-BAQ ";
	#my ($out,$stderr,$exit) = capture{system($cmd)};
	#my ($max)=$self->_parse_pileup($out);
	
}

=head2  _run_clipOver
run bamutils clipover to clip overlapping read pairs to calculate coverage peak for a gene
Inputs
=over 2
=item gene_bam bam file containing gene specific reads
=back
=cut


sub _run_clipOver {
	my($self,$gene_bam)=@_;
	my $tmp_gene_clipped=$self->options->{'tmpdir'}.'/tmp_gene_clipped';
	my $cmd="samtools view -h $gene_bam | bam clipOverlap --readName --in - --out $tmp_gene_clipped.bam";
	my ($out,$stderr,$exit)=capture{system($cmd)};
	chomp $stderr;
	if ($stderr=~m/Completed ClipOverlap Successfully/) {
		# create coordinate sorted bam ...
		Bio::DB::Bam->sort_core(0,$tmp_gene_clipped.'.bam',$tmp_gene_clipped.'_coord');
		Bio::DB::Bam->index_build($tmp_gene_clipped.'_coord.bam');
		return $tmp_gene_clipped.'_coord.bam';
	}
	else {
		$log->logcroak("ClipOverlap failed to run  OUT:$out  :ERR: $stderr EXIT:$exit");
	}
}

# parse pilpeup output
# currently no in use
sub _parse_pileup {
	my ($self,$pileup_out)=@_;	
	my $max;
	foreach my $line((split("\n", $pileup_out))) 
	{
		my ($cov)=(split "\t" , $line) [3];
		$max = $cov if $cov > $max;
	}
	$max;
}


=head2  _print_peak
print peak data to file
Inputs
=over 2
=item hash - has string peak data 
=item fh - file handler to write data
=back
=cut

sub _print_peak {
	my ($self,$hash,$fh)=@_;
	foreach my $key (sort keys %$hash) {
		if(defined($hash->{$key}))
		{	
			print $fh "$key\t".$hash->{$key}."\n";	
		}
	}
}


sub options {
shift->{'options'};
}



1;