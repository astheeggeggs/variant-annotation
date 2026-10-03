use strict; use Bio::DB::BigWig 'binMean'; use Time::HiRes qw(time);
my $f = shift;
for my $iv (["21", 33025935, 33026100], ["21", 33025935, 33040000], ["21", 46000000, 46000050]) {
  my $t = time;
  my $wig = Bio::DB::BigWig->new(-bigwig => $f);
  my @feats = $wig->features(-type => 'summary', -seq_id => $iv->[0], -start => $iv->[1], -end => $iv->[2]);
  my $tot = 0; $tot += $_->length * binMean($_->score()) for @feats;
  printf "%s:%d-%d  sum=%.6f  (%.1fs)\n", @$iv, $tot, time - $t;
}
