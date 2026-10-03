# Minimal Bio::Perl for LOFTEE. LOFTEE does `use Bio::Perl` and only calls reverse_complement(), but the
# conda BioPerl 1.7.8 packages no longer ship Bio::Perl, so without this file the LoF plugin fails to
# compile (and VEP carries on without it). The subs below are copied verbatim from BioPerl 1.6.924
# (the version the original vep105_loftee image installed via vep_install), Bio/Perl.pm, which is
# free software under the same terms as Perl itself.
package Bio::Perl;

use strict;
use Carp;
use Bio::PrimarySeq;
use base qw(Exporter);

our @EXPORT = qw(reverse_complement revcom reverse_complement_as_string revcom_as_string);
our @EXPORT_OK = @EXPORT;

sub reverse_complement {
    my ($scalar) = shift;

    my $obj;

    if( ref $scalar ) {
        if( !$scalar->isa("Bio::PrimarySeqI") ) {
            confess("Expecting a sequence object not a $scalar");
        } else {
            $obj= $scalar;
        }

    } else {

        # check this looks vaguely like DNA
        my $n = ( $scalar =~ tr/ATGCNatgcn/ATGCNatgcn/ );

        if( $n < length($scalar) * 0.85 ) {
            confess("Sequence [$scalar] is less than 85% ATGCN, which doesn't look very DNA to me");
        }

        $obj = Bio::PrimarySeq->new(-id => 'internalbioperlseq',-seq => $scalar);
    }

    return $obj->revcom();
}

sub revcom {
    return &Bio::Perl::reverse_complement(@_);
}

sub reverse_complement_as_string {
    my ($scalar) = shift;
    my $obj = &Bio::Perl::reverse_complement($scalar);
    return $obj->seq;
}

sub revcom_as_string {
    my ($scalar) = shift;
    my $obj = &Bio::Perl::reverse_complement($scalar);
    return $obj->seq;
}

1;
