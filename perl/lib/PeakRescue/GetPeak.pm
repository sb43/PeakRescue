
package PeakRescue::GetPeak; 
use PeakRescue;
our $VERSION = PeakRescue->VERSION;

use Bio::DB::Sam;
use Bio::DB::Sam::Constants;

use FindBin qw($Bin);
use List::Util qw(max min);
use Capture::Tiny qw(:all);
use Data::Dumper;
use Log::Log4perl;
Log::Log4perl->init("$Bin/../config/log4perl.gt.conf");
my $log = Log::Log4perl->get_logger(__PACKAGE__);

sub new {
	my ($class,$options)=@_;
	$log->trace('new');
	my $self={};
	bless $self, $class;
	$self->_init($options);
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
	$options->{'f'}=$options->{'o'}.'/'.$file_name;	
	$self->{'options'} = $options;
	
	return;
}

sub _get_bam_object {
	my ($bam_files,$options)=@_;
	my (%bam_objects,%bas_files);
	foreach my $file (@$bam_files) {
		my $sam = Bio::DB::Sam->new(-bam => "$options->{'d'}/$file.bam",
															-fasta =>$options->{'g'},
															-expand_flags => 1);
		$sam->max_pileup_cnt($MAX_PILEUP_DEPTH);
		$bam_objects->{$file}=$sam;
	}
	$bam_objects;
}

sub _get_max_coverage {
my ($sam_object,$g_pu)=@_;
my $new_value = 100000;
print "Chr".$g_pu->{'chr'}.'start=>'.$g_pu->{'start'}.'-end=>'.$g_pu->{'end'}."\n";
#set this value at the top
Bio::DB::Sam->max_pileup_cnt($new_value);# chanages default coverage cap from 8000 to $new_value
my ($coverage) = $sam_object->features(-type=>'coverage',-seq_id=>$g_pu->{'chr'},-start=>$g_pu->{'start'},-end=>$g_pu->{'end'}, -filter=>sub {
my $a=shift;
$a->qual > 80;

});
my @data  = $coverage->coverage;
my $max_cov=max(@data);

print "\nMAX:$max_cov\n@data\n";

}