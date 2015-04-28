package PeakRescue;
use strict;
use Const::Fast qw(const);

use base 'Exporter';
our $VERSION = '1.0.0';
our @EXPORT = qw($VERSION);
const my $LICENSE =>
"#################
# PeakRescue::GlobalTranscript version %s, Copyright (C) 2015 University of Edinburgh and Wellcome Trust Cancer Genome Project (CGP)
# This scripts comes with ABSOLUTELY NO WARRANTY
# See LICENSE for full details.
#################";

sub license {
  return sprintf $LICENSE, $VERSION;
}

1;

