package PeakRescue::RunPeakRescue; 
use PeakRescue;
our $VERSION = PeakRescue->VERSION;

use Bio::DB::Sam;
use strict;
use FindBin qw($Bin);
use List::Util qw(max min);
use File::Path qw(mkpath remove_tree);
use Capture::Tiny qw(:all);
use File::Basename;
use File::Spec;
use Data::Dumper;
use Log::Log4perl;
use Const::Fast qw(const);
Log::Log4perl->init("$Bin/../config/log4perl.gt.conf");
my $log = Log::Log4perl->get_logger(__PACKAGE__);

use PeakRescue::Base;
use PeakRescue::GlobalTranscript;
use PeakRescue::GetPeak;


# temporary paths for testing, can be set in config ini file
const my $PYTHON_PATH => "";
const my $HTSEQ_PATH => "/nfs/users/nfs_s/sb43/software/HTSeq-0.6.1p1/HTSeq/scripts";
const my $PICARD_PATH => "/software/CGP/external-apps/picard-tools-1.80/lib";

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
	$options->{'tmpdirPeak'}=$options->{'o'}.'/'.'tmp_runPeakRescue';
	
	$self->{'options'} = $options;
	
	$log->debug("Using bam file:".$options->{'bam'});
	$log->debug("Using gtf file:".$options->{'gtf'});
	$log->debug("Using genome file:".$options->{'g'});
	$log->debug("Results will be written to :".$options->{'o'});
	$log->debug("Using temp dir:".$options->{'tmpdirPeak'});
	return;
}



sub run_pipeline {
	my($self)=@_;
	$self->_run_htseq;
	$self->_run_htseq_disambiguate;
  $self->_process_sam;
  $self->_process_gtf;
  $self->_get_peak;
  $self->_runPeakrescue;
  $self->_process_output;
  
  $log->info("PeakRescue pipeline completed successfully");
  $log->info("Process log is written in peakrescue.log");
  $log->info("Results are stored in folder ".$self->options->{'o'});
  
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
	#my $gtf = $self->options->{'gtf'};	
	#my $bam = $self->options->{'bam'};
	# store output ...
	$self->options->{'htseq_sam'}=$self->options->{'tmpdirPeak'}.'/'.$self->options->{'f'}.'_htseq.sam';
	$self->options->{'htseq_count'}=$self->options->{'tmpdirPeak'}.'/'.$self->options->{'f'}.'_htseq_count.out';
	
	# requires read name sorted sam file...
	my $cmd = "samtools sort -on ".$self->options->{'bam'}." tmpsort | samtools view - | ".
		"python ".
		"$HTSEQ_PATH/count.py ".
			"--mode=union ".
			"--stranded=no ".
			"--samout=".$self->options->{'htseq_sam'}.
			" --type=exon ".
			"--idattr=gene_id ". 
			"- ".
			$self->options->{'gtf'}.
			" >".$self->options->{'htseq_count'};
	# run command
	PeakRescue::Base->_run_cmd($cmd);
	
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
	$self->options->{'disambiguated_sam'} =$self->options->{'tmpdirPeak'}.'/'.$self->options->{'f'}.'_disambiguated.sam';
	$self->options->{'disambiguated_count'}=$self->options->{'tmpdirPeak'}.'/'.$self->options->{'f'}.'_disambiguated_count.out';
	$self->options->{'multimapped_rngn'}=$self->options->{'tmpdirPeak'}.'/'.$self->options->{'f'}.'_multimapped_readname_gene_name.out';
	$self->options->{'amb_rngn'}=$self->options->{'tmpdirPeak'}.'/'.$self->options->{'f'}.'_ambiguous_readname_gene_name.out';

	my $cmd = "grep -P \"ambiguous|alignment_not_unique\" ".$self->options->{'htseq_sam'}.
		" | python ".
		"$HTSEQ_PATH/count_htseq_peakRescue.py ".
			"--mode=union ".
			"--stranded=no ".
			"--samout=".$self->options->{'disambiguated_sam'}.
			" --type=exon ".
			"--idattr=gene_id ". 
			"- ".
			$self->options->{'gtf'}." ".
			$self->options->{'multimapped_rngn'}." ".
			$self->options->{'amb_rngn'}." ".
			" >".$self->options->{'disambiguated_count'};
	# run command
	PeakRescue::Base->_run_cmd($cmd);
	
}


=head2 _process_sam
create Karyotypic sorted bam file using picard required for GATK pcoverage calculation
Inputs
=over 2
=back
=cut

sub _process_sam {
	my ($self)=@_;
  my $tmp_combined_bam=$self->options->{'tmpdirPeak'}.'/'.$self->options->{'f'}.'_combined_sorted.bam';
	$self->options->{'kayrotypic'}=$self->options->{'tmpdirPeak'}.'/'.$self->options->{'f'}.'_kayrotypic.bam';
	
  # add disambiguated reads containing additional XF:Z tags to original sam with updated 
 	my $cmd = "samtools view -H ".
 	   $self->options->{'bam'}.
 	  " | cat -  ".$self->options->{'disambiguated_sam'}. " ".$self->options->{'htseq_sam'}." | samtools view -bS -| ".
		"samtools sort -o - tmp_sort >$tmp_combined_bam";  # took 19 min to create 1.6GB file
  
    PeakRescue::Base->_run_cmd($cmd);
  
    Bio::DB::Bam->index_build($tmp_combined_bam);
    #$cmd="samtools index $tmp_combined_bam";
   	#PeakRescue::Base->_run_cmd($cmd);
   	
 # creates bam to Karyotipic sorted bam
		$cmd = "java -Xmx2G -jar $PICARD_PATH/ReorderSam.jar ".
						"I= $tmp_combined_bam ".
						"O= ".$self->options->{'kayrotypic'}.
						" REFERENCE= ".$self->options->{'g'}.
						" ALLOW_CONTIG_LENGTH_DISCORDANCE= 1";
										
		PeakRescue::Base->_run_cmd($cmd);
		Bio::DB::Bam->index_build($self->options->{'kayrotypic'});
		#$cmd="samtools index ".$self->options->{'kayrotypic'};
		#PeakRescue::Base->_run_cmd($cmd);
		
}

=head2 _process_gtf
process gtf file using GlobalTranscript package
Inputs
=over 2
=back
=cut

sub _process_gtf {
	my ($self)=@_;
	$self->options->{'u'}=1;
	my $gt=PeakRescue::GlobalTranscript->new($self->options);
	$self->options->{'unique_regions'}=$gt->{'unique_regions'};
	$self->options->{'global_transcript'}=$gt->{'global_transcript'};
	$self->options->{'geneboundaries'}=$gt->{'geneboundaries'};
	$self->options->{'global_transcript_gene_length'}=$gt->{'global_transcript_gene_length'};
	$self->options->{'unique_segment_gene_length'}=$gt->{'unique_segment_gene_length'};
	$self->options->{'non_overlapping_geneboundaries'}=$gt->{'non_overlapping_geneboundaries'};
	undef $gt;
	$self;
}

=head2 _get_peak
get peak 
Inputs
=over 2
=back
=cut

sub _get_peak {
	my ($self)=@_;
	my $getPeak_options;
	$getPeak_options->{'bed'}=$self->options->{'geneboundaries'};
  $getPeak_options->{'bam'}=$self->options->{'kayrotypic'};
  $getPeak_options->{'g'}=$self->options->{'g'};
  $getPeak_options->{'gt'}=$self->options->{'global_transcript'};
  $getPeak_options->{'o'}=$self->options->{'o'};
  $getPeak_options->{'alg'}=$self->options->{'alg'};
 	my $peak=PeakRescue::GetPeak->new($getPeak_options);
  $self->options->{'peak_file'}=$peak->{'peak_file'};
  undef $peak;
  $self;
}

=head2 _runPeakrescue
run peakrescue  
Inputs
=over 2
=back
=cut

sub _runPeakrescue {
	my ($self)=@_;
	
	#my $combinedReadNameGeneName = $self->options->{'tmpdirPeak'}.'/tmpCombinedRNGN.tab';
	#my $cmd_cat = "cat ".$self->options->{'amb_rngn'}." ".$self->options->{'multimapped_rngn'}." > $combinedReadNameGeneName";
	#PeakRescue::Base->_run_cmd($cmd_cat);
	
	# we can choose which file type to use 
	 $self->_runReadToGeneAssignment($self->options->{'multimapped_rngn'},'multimappers');
	 $self->_runReadToGeneAssignment($self->options->{'amb_rngn'},'ambiguous_unique');
	 
	 # if we would like to analyse all reads in single pass
	 #$self->_runReadToGeneAssignment($combinedReadNameGeneName,'combined');
	 
}


sub _runReadToGeneAssignment {
	my($self,$ReadNameGeneNameFile,$tag) = @_;
	my $cmd= "python $Bin/peakRescue_readToGeneAssignment.py -d ".
	$self->options->{'o'}.
	" -p ".$self->options->{'peak_file'}.
	" -m $ReadNameGeneNameFile".
	" -l ".$self->options->{'global_transcript_gene_length'}.
	" -t $tag";
	
	PeakRescue::Base->_run_cmd($cmd);
	
	#results_peakrescue_readtype_multimappers_all_genes.out
	$self->options->{$tag}=$self->options->{'o'}.'/results_peakrescue_readtype_'.$tag.'_all_genes.out';

}





=head2 _process_output
process output files to create final output file with FPKM 
Inputs
=over 2
=back
=cut

sub _process_output {
	my ($self)=@_;
	my $final_count_with_fpkm=$self->options->{'o'}.'/peakRescueFinalCount.out';
	
	# merge files in pairs using ensid as common key 
	# instead of using paste merge files by ensid key
	
	my $cmd = "$Bin/mergeFiles.pl -f1 ".$self->options->{'htseq_count'}.
	" -f2 ".$self->options->{'disambiguated_count'}.
	" -c1 1 -c2 1 -o ".$self->options->{'tmpdirPeak'}.'/tmp_htseq_count_n_disambiguated_count.out';  
	
	PeakRescue::Base->_run_cmd($cmd);
	
	$cmd = "$Bin/mergeFiles.pl -f1 ".$self->options->{'ambiguous_unique'}.
	" -f2 ".$self->options->{'multimappers'}.
	" -c1 1 -c2 1 -o ".$self->options->{'tmpdirPeak'}.'/tmp_ambiguous_unique_n_multimappers.out';
	
	PeakRescue::Base->_run_cmd($cmd);
	
	$cmd = "$Bin/mergeFiles.pl -f1 ".$self->options->{'tmpdirPeak'}.'/tmp_htseq_count_n_disambiguated_count.out'.
	" -f2 ".$self->options->{'tmpdirPeak'}.'/tmp_ambiguous_unique_n_multimappers.out'.
	" -c1 1 -c2 1 -o ".$self->options->{'tmpdirPeak'}.'/tmp_count_n_contributions.out';
	
	PeakRescue::Base->_run_cmd($cmd);
	# add length
	$cmd = "$Bin/mergeFiles.pl -f1 ".$self->options->{'tmpdirPeak'}.'/tmp_count_n_contributions.out'.
	" -f2 ".$self->options->{'global_transcript_gene_length'}.
	" -c1 1 -c2 1 -o ".$self->options->{'tmpdirPeak'}.'/tmp_count_n_contributions_n_length.out';
	
	PeakRescue::Base->_run_cmd($cmd);
	
  #gene   uc  mm_c total    gene    uc_d    mm_tr   amb_tr  gene    amb_p   gene    mm_p    ensid   gene_len_gt 
  #0      1      2   3      4        5       6       7       8       9       10      11      12     13        
  my $total_read_count;
  my $all_data;
  open (my $fh_data, '<', $self->options->{'tmpdirPeak'}.'/tmp_count_n_contributions_n_length.out');
  while (<$fh_data>) {
    chomp;
    my($gene, $uc, $mm_c, $uc_d, $amb_tr, $amb_p, $mm_p, $length)=(split "\t", $_)[0,1,2,5,7,9,11,13];
    my $final_uc=($uc+$uc_d+$amb_p+$mm_p);
    my $all_nc=($mm_c+$uc_d+$amb_tr);
    $total_read_count+=$final_uc;
    my $line;
    push(@$line,$uc,$mm_c,$all_nc,$uc_d,$amb_p,$mm_p,$final_uc,$length);
    $all_data->{$gene}=$line;
  }
  open(my $fh_final,'>', $final_count_with_fpkm);
	#fpkm loop 
	print $fh_final "Gene\tuniqueCount\tmultiMappedCount\ttotalNonUniqueCount\tuniqueDisambiguatedCount\tambiguousProportion\tmultimappedProportion\ttotalCount\tglobalTranscriptLength\tfpkmTotalCount\n";
  foreach my $gene (keys %$all_data) {
    my $fpkm = ((10**9) *  @{$all_data->{$gene}}[6] ) / ($total_read_count * @{$all_data->{$gene}}[7] );
    print $fh_final "$gene\t".join("\t",@{$all_data->{$gene}})."\t".sprintf("%.3f",$fpkm)."\n";
  }
  PeakRescue::Base->cleanup_dir($self->options->{'tmpdirPeak'});
  
  
}



sub options {
shift->{'options'};
}



1;