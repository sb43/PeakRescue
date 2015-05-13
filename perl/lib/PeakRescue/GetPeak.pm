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
use Try::Tiny qw(try catch finally);
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
	return;
}

sub _get_bam_object {
	my ($self,$bam_file)=@_;
		my $bam = Bio::DB::Sam->new(-bam => $bam_file,
															-fasta => $self->options->{'g'},
															-expand_flags => 1,
															-split_splices => 1
															);
		$bam->max_pileup_cnt($MAX_PILEUP_DEPTH); # chanages default coverage cap from 8000 to $new_value
		return $bam;
}

sub _do_max_peak_calculation { 
	my ($self)=@_;
	my($bam_object)=$self->_get_bam_object($self->options->{'bam'});
	my($tabix)=$self->_get_tabix_object();
	open(my $peak_fh, '>', $self->options->{'o'}.'/'.$self->options->{'f'}.'_peak.txt');
	$log->logcroak("Unable to crea tabix object") if (!$tabix);
	# get chromosome names
	my @chrnames=$tabix->getnames;
	###
	# Querying is ALWAYS half open regardless of underlying file type [ i.e zero start and 1 end -- same as bed file  ]
	###
	foreach my $chr (@chrnames) {
	my ($peak_data, $counter);
	my $res = $tabix->query($chr);
		while(my $record = $tabix->read($res)){		
			chomp;
			my ($chr,$start,$end,$gene)=(split "\t", $record) [0,1,2,3];
			$counter++;
			if($counter % 500 == 0) {
				$log->debug("Calculated peak for chromosome: $chr ==>".$counter.' genes');
				$self->print_peak($peak_data,$peak_fh);
				$peak_data=();
			}
				
			#my($gene_bam_file)=$self->_get_gene_bam_object2($bam_object,"$chr:$start-$end",$gene);
			
			my($max_peak)=$self->_get_gene_bam_object2($bam_object,$chr,$start,$end,$gene);
			#my($gene_bam_object) = $self->_get_gene_bam_object($bam_object,$chr,$start,$end,$gene);
			#my $segment=$gene_bam_object->segment( -seq_id => $chr, -start => $start, -end => $end);
			
			#my ($coverage) = $segment->features('coverage');
			#my ($coverage) = $gene_bam_object->features(-type => 'coverage', -seq_id => $chr, -start => $start, -end => $end);
			
			#store temporary gene interval 
			#my $tmp_interval_bed=$self->options->{'tmpdir'}.'/tmp_interval.bed';
			#open (my $tmp_int_fh, '>',$tmp_interval_bed); 
			#print $tmp_int_fh "$chr\t$start\t$end\n";
			#$self->_get_bed_coverage($gene_bam_file,$tmp_interval_bed);
			#	exit;
			#print "$gene: $counter----> $max_peak \n";
			$peak_data->{$gene}=$max_peak;
		  #undef $gene_bam_object;
		}
	$self->print_peak($peak_data,$peak_fh);
	$log->debug(">>>>>>>>> Completed peak calculation for chromosome: $chr ===>".$counter.' gene >>>>>>>>>');
}
	close($peak_fh);
	
}

sub _get_tabix_object {
	my ($self)=@_;
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


sub _get_gene_bam_object {
	my ($self,$bam_object,$chr,$start,$end,$gene)=@_;
	my $tmp_gene_file=$self->options->{'tmpdir'}.'/tmp_gene.bam';
	my $bam = Bio::DB::Bam->open($tmp_gene_file,'w');
	$bam->header_write($bam_object->header);
	my @alignments = $bam_object->get_features_by_location( -seq_id => $chr,
																									 -start  => $start,
																									 -end    => $end);
																									 																							 
	 for my $align (@alignments) {
			my $str = $align->aux;
			if ($str =~m/XF:Z:$gene/) {
			  	# alignments filtered  by given gene 
			  $bam->write1($align->{'align'});
			}
	 }
	# this is required to put EOF marker before creating bam object
	undef $bam;
	#system("samtools index $tmp_gene_file");
	Bio::DB::Bam->index_build($tmp_gene_file);
	#my ($gene_bam)=$self->_get_bam_object($tmp_gene_file);
	return $tmp_gene_file;
}

=head2  _get_gene_bam_object2
fetch reads
Inputs
=over 2
=item sam_object - Bio::DB sam object
=item region - chr:start-stop format region info to get reads
=item Reads_FH - temp file handler to store reads
=back
=cut

sub _get_gene_bam_object2 {
my ($self,$bam_object,$chr,$start,$end,$gene)=@_;
	my $tmp_gene_file=$self->options->{'tmpdir'}.'/tmp_gene.bam';
	my $bam = Bio::DB::Bam->open($tmp_gene_file,'w');
	$bam->header_write($bam_object->header);
	my $tmp_gene_sorted=$self->options->{'tmpdir'}.'/tmp_gene_name_sorted';
	my $read_flag=undef;
	$bam_object->fetch("$chr:$start-$end", sub {
		my $a = shift;
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
	if (!$read_flag ) {
		#print "No reads for :  $gene\n";
		return "0";
	}
	
	#system("samtools index $tmp_gene_file");
	Bio::DB::Bam->index_build($tmp_gene_file);
	# need samtools 1.1 or above which takes care of overlapping read pairs...
	# sort bam file
	##Bio::DB::Bam->sort_core(1,$tmp_gene_file,$tmp_gene_sorted);
	#my($clipped_bam)=$self->_run_clipOver("$tmp_gene_sorted.bam");
	my ($gene_bam)=$self->_get_bam_object($tmp_gene_file);
  my ($coverage) = $gene_bam->features(-type => 'coverage', -seq_id => $chr, -start => $start, -end => $end);
  #option2 samtools mpipeup -- very slow/ GATK very slow...
	#my $cmd= "~/software/samtools-1.2/samtools mpileup $tmp_gene_file -d $MAX_PILEUP_DEPTH -A -f ".$self->options->{'g'}. " -r $region --no-BAQ ";
	#my ($out,$stderr,$exit) = capture{system($cmd)};
	#my ($max)=$self->_parse_pileup($out);
	
	return max($coverage->coverage);
}

sub _run_clipOver {
	my($self,$gene_bam)=@_;
	my $tmp_gene_clipped=$self->options->{'tmpdir'}.'/tmp_gene_clipped';
	my $cmd="samtools view -h $gene_bam | bam clipOverlap --readName --noPhoneHome --in - --out $tmp_gene_clipped.bam";
	
	my ($out,$stderr,$exit)=capture{system($cmd)};
	chomp $stderr;
	if ($stderr=~m/Completed ClipOverlap Successfully/) {
		Bio::DB::Bam->sort_core(0,"$tmp_gene_clipped.bam","$tmp_gene_clipped\_name");
		Bio::DB::Bam->index_build("$tmp_gene_clipped\_name.bam");
		return "$tmp_gene_clipped\_name.bam";
	}
	else {
		$log->logcroak("ClipOverlap fialed to run OUT:$out  :ERR: $stderr EXIT:$exit");
	}
	
}


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

sub print_peak {
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