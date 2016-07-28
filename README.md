# chipfiring
This is a collection of files related to the master thesis "Chip Firing, Riemann-Roch, and Brill-Noether Theory for Graphs", by Alejandro Vargas. The contents are as follows:

The main file is ChipFiring.pm, a perl module implementing: Dhar's burning algorithm, Iterated Dhar's burning algorithm for calculating reduced divisors, a function that checks whether a divisor is rank-1, and a function that exhaustively finds all the rank-1 divisors. The module has several dependencies which may be installed using cpan.

A folder of examples with several programs that the thesis makes reference to.

A folder of extra files with a data set and a copy of the thesis. In it there is a description of the algorithm used to find all the rank-1 divisors.

We illustrate the main methods of ChipFiring.pm with a simple program. Consider a wedge of two triangles (see example 1.16). It has five vertices numbered from 0 to 4 and edges "0,1 0,2 0,3 0,4 1,2 3,4". This graph has a single rank-1 divisor of degree 2, which is found by the all rank-1 function. The all rank-1 function searchs for rank-1 divisors of degree floor(g/2) + 1, where g is the genus of the graph. This is done in order to test a conjectured Brill-Noether theory for graphs. We reduce this divisor with respect to v4.



```perl
use strict;
use warnings;
use v5.20;
use ChipFiring;

my $graph = ChipFiring->new();
$graph->from_links("0,1 0,2 0,3 0,4 1,2 3,4");
say "The adjacency matrix is", "\n", $graph->matrix_string();
$graph->all_rank1();
my @divisors = $graph->rank1_string();

say "The rank-1 divisors are", "\n", join "\n", @divisors;
say "The divisor 2v0 has rank ", $graph->rank("2v0");
say "The reduction of the first rank-1 divisor with respect to v4 is ", $graph->reduce(4, $divisors[0]);
```



Running the program yields the following output:
```
The adjacency matrix is
0	1	1	1	1
1	0	1	0	0
1	1	0	0	0
1	0	0	0	1
1	0	0	1	0
The rank-1 divisors are
2v0
The divisor 2v0 has rank 1
The reduction of the first rank-1 divisor with respect to v4 is 1v3 1v4
```


