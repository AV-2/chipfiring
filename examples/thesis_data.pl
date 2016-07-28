use strict;
use warnings;
use v5.20;
use ChipFiring;

my $g = $ARGV[0];
my $n = $ARGV[1];
my $tries = $ARGV[2];
my $tries2 = $ARGV[3];
my $jump = $ARGV[4];

ChipFiring->use_database('/home/av/temp/thesisdata.sqlite', "g$g" . "v$n");
$ChipFiring::MAX_ATTEMPTS = 20;

for (1..500) {
say $n;
for (1..$tries2) {
    
    my $string = '';
    for (1..$tries) {
        my $graph = ChipFiring->new();
        $graph->random_graph($n, $g);
        $graph->all_rank1();
        $string .=  $graph->rank1_db_string() . ',';
    }
    #remove the last comma
    chop $string; 
    ChipFiring->commit_database('rank1', $string);
}
$n += $jump;
ChipFiring->use_table( "g$g" . "v$n");

}

