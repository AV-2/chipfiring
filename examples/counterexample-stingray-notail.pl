use strict;
use warnings;
use v5.20;
use ChipFiring;


my $g = 6;
my $n = 8;

my $matrix = <<'MATRIX';
0,1
0,6
0,6
1,5
1,3
2,1
2,1
3,4
3,4
4,5
4,6
5,6
MATRIX



my $graph = ChipFiring->new();
$graph->from_links($matrix);

my $graph2 = ChipFiring->new(signature => 'sig');
$graph2->from_links($graph->subdivide(1));
say $graph2->matrix_string();
$graph2->all_rank1();
my $graph_string =  $graph2->rank1_db_string() . ',';

#Let us load the divisors into a hash 
my %divisors;
$divisors{$_} = 1 foreach @{$graph2->rank1_divisors()};

#We get divisors with points inside edges. We know that in our algorithm points inside edges have an index above n
#say $_ for (@{$graph2->rank1_divisors()}); 
my @divisors = grep {/(\d+)$/; $1 >= $n} @{$graph2->rank1_divisors()};

#say $_, " ", join '-', $graph2->neighborhood($_) foreach (@divisors);
      
foreach my $div (@divisors) {
    my $rank1_approx = count_rank_approximations($graph2, $div, \%divisors); 
    if ($rank1_approx == 0) {
        say "We found a problem";
        say $div, "\n";
        
        say "The rank 1 divisors are";
        say join "\n", @{$graph2->rank1_divisors()}, "\n";
        
        #get all the equivalent divs
        my @equivalent_divs = $graph2->all_equivalent($div);
        say "The equivalent effective divisors are";
        say join "\n",@equivalent_divs, "\n";
        
        #count
        foreach my $div (@equivalent_divs) {
            say "We now approximate $div";
            my @approx = $graph2->div_approximate($div, $n);

            #make them 0 reduced
            say "For ", $approx[$_], " it has a v0-reduced divisor ", $graph2->reduce(0, $approx[$_]), " with rank ", (defined($divisors{$graph2->reduce(0, $approx[$_])})) ? 1:0 for (0..$#approx);
        }
        say "\n";
        say $graph2->links_string();
        say $graph2->matrix_string();
    }

    
}

sub count_rank_approximations {
    my $graph = shift;
    my $div = shift;
    my %divisors = %{shift()};
    
    my $rank1_approx = 0;
    
    #running the get all equivalent divisors routine is expensive
    #So we first approximate the divisor itself, before approximating all its equivalent ones.
    my @approx = $graph->div_approximate($div, $n);
        
    #make them 0 reduced
    $approx[$_] = $graph->reduce(0, $approx[$_]) for (0..$#approx);
    $rank1_approx += defined($divisors{$_}) foreach (@approx);
    
    return $rank1_approx if ($rank1_approx);
    
    #get all the equivalent divs
    my @equivalent_divs = $graph->all_equivalent($div);
    
    #count
    foreach my $div (@equivalent_divs) {
        my @approx = $graph->div_approximate($div, $n);
        
        #make them 0 reduced
        $approx[$_] = $graph->reduce(0, $approx[$_]) for (0..$#approx);
        
        $rank1_approx += defined($divisors{$_}) foreach (@approx);
        last if ($rank1_approx);
    }
    
    return $rank1_approx;
}

sub get_equivalent_{
    my $graph = shift;
    my $div = shift;
    my %divisors = %{shift()};
    
    my @approx = $graph->div_approximate($div, $n);
    $approx[$_] = $graph->reduce(0, $approx[$_]) for (0..$#approx);
    my $rank1_approx = 0;
    $rank1_approx += defined($divisors{$_}) foreach (@approx);
    return $rank1_approx;
}