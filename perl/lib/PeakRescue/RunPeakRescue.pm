package PeakRescue::RunPeakRescue; 
use PeakRescue;
our $VERSION = PeakRescue->VERSION;

use strict;
use FindBin qw($Bin);
use List::Util qw(max min);
use File::Path qw(mkpath remove_tree);
use Capture::Tiny qw(:all);
use File::Basename;
use File::Spec;
use Data::Dumper;
use Log::Log4perl;
Log::Log4perl->init("$Bin/../config/log4perl.gt.conf");
my $log = Log::Log4perl->get_logger(__PACKAGE__);

use PeakRescue::Base;

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
	mkpath($options->{'o'}.'/'.'tmp_runPeakRescue');
	$options->{'tmpdir'}=$options->{'o'}.'/'.'tmp_runPeakRescue';
	
	$self->{'options'} = $options;
	
	$log->debug("Using bam file:".$options->{'bam'});
	$log->debug("Using gtf file:".$options->{'gtf'});
	$log->debug("Using genome file:".$options->{'g'});
	$log->debug("Results will be written to :".$options->{'o'});
	$log->debug("Using temp dir:".$options->{'tmpdir'});
	return;
}



sub run_peakrescue {
	my($self)=@_;
	$self->_run_htseq;
	$self->_run_htseq_disambiguate;
  $self->process_sam;
  
  # run rest of the steps
  
  $log->info("PeakRescue pipeline completed successfully");
  $log->info("Process log is written in peakrescue.log");
  $log->info("Results are stored in ".$self->options->{'o'});
  
  return 1;
}


=head2 _run_htseq
run htseq modified version which outputs sam file with XF:Z flag and counts data containing multimapped reads ?? to be confirmed
Inputs
=over 2
=item gtf - GTF file
=item bam - bam file from splice aware aligner 
=back
=cut


sub _run_htseq {
	my($self)=@_;
	my $gtf = $self->options->{'gtf'};	
	my $bam = $self->options->{'bam'};
	my $htseq_sam=$self->options->{'tmpdir'}.'/'.$self->options->{'f'}.'_htseq.sam';
	my $htseq_count=$self->options->{'tmpdir'}.'/'.$self->options->{'f'}.'_htseq_count.out';
	
	# requires read name sorted sam file...
	my $cmd = "samtools sort -on $bam tmpsort | samtools view - | ".
		"python ".
		"/nfs/users/nfs_s/sb43/software/HTSeq-0.6.1p1/HTSeq/scripts/count_0.5.3p3.py ".
			"--mode=union ".
			"--stranded=no ".
			"--samout=$htseq_sam ".
			"--type=exon ".
			"--idattr=gene_id ". 
			"- ".
			"$gtf ".
			">$htseq_count";
	  
	my ($out,$stderr,$exit)=capture{system($cmd)};
	if($exit) {
			$log->logcroak("HTSeq count step failed with status <<<<<< OUT:$out  :ERR: $stderr EXIT:$exit");
	}
	else {
		$log->debug("Step1 <<<<<<< HTSeq completed successfully");
	}
	
	$self->options->{'htseq_sam'}=$htseq_sam;
	$self->options->{'htseq_count'}=$htseq_count;
}

=head2 _run_htseq_disambiguate
run htseq modified version to outputs sam file with additional XF:Z flag for disambiguated reads  
and counts data containing unique disambiguated reads , ambiguous reads and multimapped reads 
Inputs
=over 2
=item gtf - GTF file
=item htseq_sam - sam file from previous step
=back
=cut

sub _run_htseq_disambiguate {
	my($self)=@_;
	my $gtf = $self->options->{'gtf'};	
	my $htseq_sam=$self->options->{'htseq_sam'};
	my $disambiguated_sam=$self->options->{'tmpdir'}.'/'.$self->options->{'f'}.'_disambiguated.sam';
	my $disambiguated_count=$self->options->{'tmpdir'}.'/'.$self->options->{'f'}.'_disambiguated_count.out';
	my $multimapped_rngn=$self->options->{'tmpdir'}.'/'.$self->options->{'f'}.'_multimapped_readname_gene_name.out';
	my $amb_rngn=$self->options->{'tmpdir'}.'/'.$self->options->{'f'}.'_ambiguous_readname_gene_name.out';

	# requires read name sorted sam file...
	my $cmd = "grep -P \"ambiguous|alignment_not_unique\" $htseq_sam | ".
		"python ".
		"/nfs/users/nfs_s/sb43/software/HTSeq-0.6.1p1/HTSeq/scripts/count_htseq_v31.py ".
			"--mode=union ".
			"--stranded=no ".
			"--samout=$disambiguated_sam ".
			"--type=exon ".
			"--idattr=gene_id ". 
			"- ".
			"$gtf ".
			"$multimapped_rngn ".
			"$amb_rngn ".
			">$disambiguated_count";
	
	my ($out,$stderr,$exit)=capture{system($cmd)};
	if($exit) {
			$log->logcroak("HTSeq disambiguation step failed with status <<<<<< OUT:$out  :ERR: $stderr EXIT:$exit");
	}
	else {
		$log->debug("Step2 <<<<<<< HTSeq disambiguation step completed successfully");
	}
	
	$self->options->{'disambiguated_sam'}=$disambiguated_sam;
	$self->options->{'disambiguated_count'}=$disambiguated_count;
	$self->options->{'multimapped_rngn'}=$multimapped_rngn;
	$self->options->{'amb_rngn'}=$amb_rngn;
}

sub process_sam {

 # creates sam to Karyotipic sorted bam

}



sub combine_rngn {
# combine multimapped_rngn and amb_rngn  OR output single file from modified htseq whichever is quicker


}


sub process_gtf {

 # call processGTF.pl

}

sub get_peak {

 # call

}


sub runPeakrescue {

# python peakRescue_readToGeneAssignment.py -p peak.tab -r readNameGeneName.tab 
# output peakContibutions.out
}


sub process_output {

# mergeFiles.py -bed geneboundaries.bed -f1 htseq_count.out -f2 disambiguated_count.out -f3 peakContibutions.out  -o peakRescueFinalCount.out

# output PeakRescueFinalResults.out

}


sub options {
shift->{'options'};
}



1;