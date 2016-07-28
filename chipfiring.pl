use strict;
use warnings;
use v5.20;
use ChipFiring;
use Data::Dumper;

my $graph;
print 'chipfiring: ';

while (my $input = <STDIN>) {
    chomp($input);
    
    
    given ($input) {
        when (/matrix-file/) {$graph = ChipFiring->new()->from_matrix_file((split /\s/, $input)[1]);}
        when (/links-file\s+(.*)/) {$graph = ChipFiring->new()->from_links_file($1);}
        when ('links') {say $graph->links_string()}
        when (/reduce\s+(\d+)\s+(.*)/) {say join "\n", @{$graph->reduce_trace($1, $2)}}
        default { say $input;}
    }
    
    print "\nchipfiring: ";
}