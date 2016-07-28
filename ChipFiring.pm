package ChipFiring;
use v5.20;
use strict;
use warnings;

use Math::MatrixReal;
use Data::Dumper;
use Math::Prime::Util qw/factor_exp gcd/;
use Carp;
use Math::Round qw/round/;
use POSIX qw/ceil/;
use DBI;
use Time::HiRes qw(gettimeofday tv_interval);

=head1 Package Parameters
-\$ChipFiring::MAX\_ATTEMPTS = 13;

=head1 Attributes

A ChipFiring object has the following attributes which can be get as OBJ->attribute(); It is usually unsafe to set directly this attributes.

B<v, e, g>: If the ChipFiring object is holding a graph $G$, then $G$ has $n$ vertices, $m$ edges and genus $g$.

B<matrix>: The $n$ by $n$ adjacency matrix $A$ for the graph $G$ with vertices $v_0, v_1, v_2\ldots$

B<matrix\_string>: A string representation of the matrix.

B<lap>: The Laplacian matrix $L = D - A$ of $G$, where $D$ is a diagonal matrix encoding the degrees of the vertices and $A$ is the adjacency matrix.

B<support>: The neighborhood of vertex $v_i$.

B<support\_lap>: The indices of non-zero entries for the i-th row of the laplacian. Indexing begins at 0.

B<vertices>: An array with the indices of the vertices of the graph.

B<det>: The amount of Spanning Trees in the graph, also corresponding to the order of the Jacobian.

B<Li>: Let $L$ be the laplacian of $G$. One can get a generalized inverse of $L$, that is matrix $M$ such that $LML = L$, by deleting the $i$-th row and $i$-th column from $L$, then inverting the resulting matrix, and finally inserting a column and a row full of zeroes at the $i$-th position. This matrix is usually denoted by $L_i$. For energy calcluations we fix a vertex $q$, and instead of working with $L_q$ we work with det $\times L_q$, where det is the number of spanning trees, and the arithmetic is done modulo det.

B<q>: The current distinguished vertex for energy calculations.

B<ginverses>: Hash containing other generalized inverses. In general every time an $L_i$ gets calculated, it is thrown into this hash using $i$ as a key. Might contain other exotic generalized inverses, maybe the penrose generalized inverse for example. 

B<Li\_basis>: A basis to be used with an $L_i$ matrix. It is particularly important the mention of the $L_i$ matrix they are associated to because the code omits the $i$-th coordinate of the basis, as it vanishes under operations with the matrix $L_i$. Recovering this coordinate is pretty simple, as the degree of the basis element must be zero. 

B<q\_basis>: The current $L_i$-basis that corresponds to the current $q$ value.
=cut

use Class::Tiny qw/
    n m g c
    matrix
    matrix_string
    edges
    neighbors
    links_string
    degrees
    lap
    support
    support_lap
    vertices
    vertices_sorted
    det
    Li
    q
    ginverse
    rank1
    rank1_string
    signature
    rank1_db_string
    prime_string
    rank1_divisors
    empty_divisor
    /, {signature => ''};
    
our $MAX_ATTEMPTS = 10;
our $n;
our @mat;
our @lap;
our @support;
our @support_lap;
our @vertices;
our @vertices_sorted;
our @empty_divisor;
our @degrees;
our $c;
our $dbh;
our $table_name;

#=====================Graph processing===================

sub _process_graph {
=head1 _process_graph
    
    Handles populating the different properties of a graph. 
    
=cut
    my $self = shift();
    
    #Alias Mat, figure out v, create vertex array
    @mat = @{$self->{matrix}};
    $self->{n} = scalar @mat  if (!$self->{n});
    $self->{vertices} = [(0 .. $self->{n}-1)];
    $self->{empty_divisor} = [(0) x $self->{n}];  
    $self->{support} = [map {my @row = @{$_}; [grep {$row[$_]} @{$self->{vertices}}]} @mat];
    if (!$self->{degrees}) {
        $self->{degrees} = [];
        foreach my $i (@{$self->{vertices}}) {
            $self->{degrees}[$i] += $mat[$i][$_] foreach (@{$self->{support}[$i]});
        }
    }
    
    #calculate the edges and the genus
    $self->{m} = 0;
    $self->{m} += $_ foreach (@{$self->{degrees}});
    $self->{m} /= 2; 
    $self->{g} = $self->{m} - $self->{n} + 1;
    $self->{c} = $self->{g}/2 + 1;

    #sort the degrees 
    $self->{vertices_sorted} = [sort { -1 * ($self->{degrees}[$a] <=> $self->{degrees}[$b]) }  @{$self->{vertices}}]; # use numeric comparison
    
    #Calculate laplacian
    for my $i (@{$self->{vertices}}) {
        $self->{lap}[$i] = [map {-$_} @{$mat[$i]}];
        $self->{lap}[$i][$i] = $self->{degrees}[$i]; 
    }
    $self->{support_lap} = [map {my @row = @{$_}; [grep {$row[$_]} @{$self->{vertices}}]} @{$self->{lap}}];
    
    
    $self->{matrix_string} = _matrix2string(@{$self->{'matrix'}});
    
    
    
    #Make a link representation
    $self->_matrix2links();
    $self->{links_string} = _links2string(@{$self->{edges}});
}

sub _matrix2links {
=head1 _matrix2links
    
    Handles making a (vertex, vertex) representation of the edges of the graph
    
=cut
    my $self = shift();
    $self->_context();

    my @edges;
    my @neighbors;
    push @neighbors, [] for (1..$n);
    
    #iterate over the adjacency matrix, adding links
    my $i = 0;
    for my $vertex (@vertices) {
        for my $j ($vertex+1..$n-1) {
            if ($mat[$vertex][$j] > 0) {
                my @links;
                push @links, [$i++, $vertex, $j, 0, 0] for (1..$mat[$vertex][$j]);

                push @edges, @links;
                push @{$neighbors[$vertex]}, @links;
                push @{$neighbors[$j]}, @links;
            }
        }
    }
    
    $self->{edges} = \@edges;
    $self->{neighbors} = \@neighbors;

    _clean_context();
}

sub from_matrix {
    my $self = shift;
    $self->{'matrix'} = shift;
    
    #now call the matrix processing
    $self->_process_graph();
}

sub from_links {
    my $self = shift;
    my @links = @_;
    
    #convert any line break into spaces
    $links[$_] =~ s/\n/ /g foreach (0..$#links);
    
    $self->from_matrix(_links2matrix(@links));
    
    return $self;
}

sub from_links_file {
    my $self = shift;
    my $file = shift;
    
    open my $fh, '<', $file or croak("The provided file is inexistent");
    
    my $link_string;
    while (<$fh>) {
        chomp($_);
        $link_string .= $_ ;
    }
    $self->from_links($link_string);
    
    return $self;
}

sub from_matrix_file {
    my $self = shift;
    my $file = shift;
    
    open my $fh, '<', $file or croak("The provided file is inexistent");
    
    my $matrix_string;
    while (<$fh>) {
        chomp($_);
        $matrix_string .= $_ . "\n";
    }
    $self->from_matrix_string($matrix_string);
    
    return $self;
}

sub from_matrix_string {
    my $self = shift;
    
    $self->from_matrix(_string2matrix($_[0]));
}

#================================Graph Generation================================

sub random_graph {
=pod
B<Random Graph>: Produce a random graph. There are several random graph methods available: 
    
    !!to be implemented
=cut

    my $self = shift;
    my $n = $self->{n} = shift;
    my $g = $self->{g} = shift;
    my $method = shift || 1;
    
    #We enclose it in an infinite loop until we get connectivity, or we hit the value MAX ATTEMPTS 
    my $i = 0;
    while (1) {
        randomgraph12($self, $n, $g, 1) if ($method == 1);

        last if ($self->is_connected());
        croak("The requested parameters made us reach the number of Max Attempts before finding a connected graph. Probably v is very high in relation to g, consider raising the variable ChipFiring::MAX_ATTEMPTS if you want me to try harder to find a graph, or modifying the random graph generation method. Default value ten attempts. Attempted $i times.") if (++$i >= $MAX_ATTEMPTS);
    }

    return $self;
}

sub randomgraph12 {
    #Hey, who would have guessed writing a \emf{proper} random graph subroutine would be such a pain.
    #this is method 1 and 2, depending on the connectivity of the graph
    my $self = shift();
    my $n  = shift;
    my $g =  shift;
    my $method = shift;
    my $m = $self->{m} = $g + $n - 1;
    my (@mat, $r1, $r2, $i1, $i2);
    
    my $min_connectivity = $method; #You would need to change this calculation if you wish to mess with the order of matrix generating methods.
    
    my @to_connect = (0 .. $n-1);
    my @degrees = (0) x $n;
    push @mat, [(0) x $n] for (1..$n);

    #pick a random pair and connect them    
    for (1..$m) {
        if (@to_connect > 1) {
            my $n = @to_connect;
            
            #pick a random number and insert the edge
            $r1 = int(rand($n));
            $r2 = $r1;
            $r2 = int(rand($n)) while ($r1 == $r2);
            
            #order them first the biggest one, so the array splice later doesn't get fucked up
            ($r1, $r2) = ($r2, $r1) if ($r1 < $r2);
            
            $i1 = $to_connect[$r1];
            $i2 = $to_connect[$r2];
            
            $mat[$i1][$i2] += 1;
            $mat[$i2][$i1] += 1;
            $degrees[$i1]++; $degrees[$i2]++;
            
            splice @to_connect, $r1, 1  if ($degrees[$i1] > $method);
            splice @to_connect, $r2, 1  if ($degrees[$i2] > $method);

        } elsif (@to_connect == 1) {
            $r1 = $to_connect[0];
            $r2 = $r1;
            $r2 = int(rand($n)) while ($r1 == $r2);
            
            #add them to the matrix 
            $mat[$r1][$r2] += 1;
            $mat[$r2][$r1] += 1;
            
            #increase the degree
            $degrees[$r1]++; $degrees[$r2]++;
            
            @to_connect = () if ($degrees[$r1] > $method);
        }
        
        else {
            $r1 = int(rand($n));
            $r2 = $r1;
            $r2 = int(rand($n)) while ($r1 == $r2);
            
            #add them to the matrix 
            $mat[$r1][$r2] += 1;
            $mat[$r2][$r1] += 1;
            
            #increase the degree
            $degrees[$r1]++; $degrees[$r2]++;
        }
        
        #Some tests
        #say _matrix2string(@mat);
        #say join " ", @to_connect;
        #say join " ", @degrees;
        #say "========\n";
    }
    
    $self->{'matrix'} = \@mat;
    $self->{'degrees'} = \@degrees;
    
    #now call the matrix processing
    $self->_process_graph();
    
    return $self;
}

sub subdivide {
=head1 _divide_edge
    
    Subdivide all the edges in the graph by the specified number, return the link specification
    
=cut

    my $self = shift;
    my $a = shift;
    
    my $string;
    my $i = $self->{n};
    
    foreach my $edge (@{$self->{edges}}) {
        $string .= "$edge->[1],$i\n";
        
        if ($a > 1) {
            foreach (2..$a) {
                $i++;
                $string .= $i - 1 . ",$i\n"
            }
        }
        
        $string .= "$i,$edge->[2]\n";
        $i++;
    }
    
    return $string;
}

sub subdivide_edge {
=head1 _divide_edge
    
    Giving the id of an edge, insert the specified number of vertices inbetween the edge, return the resulting matrix
    !!!Not implemented
    
=cut
}

#======================Chip Firing Methods=======================
sub all_rank1 {
    my $self = shift;
    my ($div_str, $rank1, $rank2) = ('', 0, 0);
    $self->_context();
    
    my $divisors =  _all_divisors_iterator1($c - 1);
    $self->{rank1} = $divisors->{rank1}; 
    $self->{rank1_divisors} = [map {$_->[0]} @{$divisors->{rank1}}]; 
    $self->{rank1_string} = join "\n",
        map {
            $rank1++ if ($_->[1] == 1);
            $rank2++ if ($_->[1] > 1);
            $_->[0]
            } @{$divisors->{rank1}};
    $self->{time_interval} = $divisors->{time_interval}; 
    $self->{reduced_divisors} = $divisors->{reduced_divisors};
    
    #Get some extra info 
    my $link_string = $self->{links_string}; $link_string =~ s/\n//g;
    $self->{det} = round(abs(Math::MatrixReal->new_from_cols(\@lap)->minor(1,1)->det())) unless ($self->{det}); 
    $self->{prime_string} = join '*', map {join '^', @{$_}} factor_exp($self->{det});
    
    #Put together a database string
    #We call signature with the getter so we can get the default value
    $self->{rank1_db_string} =  sprintf("('%s', '%s', %s, %s, %s, %s, %s, %s, '%s', %s, '%s', '%s')", $link_string, $self->{rank1_string}, $rank1, $rank2, $rank1 + $rank2, $self->{reduced_divisors}, $self->{g}, $self->{n}, $self->signature(), $self->{det}, $self->{prime_string}, $self->{time_interval}); 
    
    _clean_context();
    
    return $divisors;
}

sub all_equivalent {
    my $self = shift;
    my $div_str = shift;
    $self->_context();
    
    $div_str = array2div(_reduce(0, [_div2array($div_str)])  );
    return _all_divisors_iterator2($c, \&ChipFiring::_equivalent_to_reduced, $div_str, 0); 
    
    _clean_context();
    
    #return $divisors;
}

sub _equivalent_to_reduced {
=head1 _equivalent_to_reduced

Compares the first argument, divisor in vector form, to the second argument, a reduced divisor with respect to q, the third argument in string representation. Returns a string representation of the first divisor if they are equivalent, otherwise returns false.

Call with context on.
=cut
    my @div = @{$_[0]};
    my $div = $_[1];
    my $q = $_[2];
    
    my $div_string = array2div(@div);
    my $candidate_div  = array2div(_reduce($q, \@div));
    my $rank;

    return $div_string if ($div eq $candidate_div);
    return 0;
}

sub rank {
    my $self = shift;
    $self->_context();
    
    my $rank = _rank([_div2array(shift())])->[0];
    
    _clean_context();
    
    return $rank;
}

sub reduce {
    my $self = shift;
    my $q = shift;
    $self->_context();
    
    my $reduce = array2div(_reduce($q, [_div2array(shift())]));
    
    _clean_context();
    
    return $reduce;
}

sub reduce_trace {
    my $self = shift;
    my $q = shift;
    $self->_context();
    
    my $steps = (  _reduce_footprint($q, [_div2array(shift())], 1)   )[2];
    
    _clean_context();
    
    return  $steps;
}

sub _rank {
=head1 rank

Given an effective divisor, return 0 if rank 0, 1 if rank 1, 2 if rank above 1.
Work on rank determining sets could make this more efficient, but we would have to write code detecting the topology of the graph.

=cut
    my @div = @{shift()};
    my @footprint = map {$_ > 0} @div;
    my $rank = 2;
    my $q = 0;
    
    foreach my $i (@vertices_sorted) {
        next if ($footprint[$i]);
        @_ = _reduce_footprint_stop($i, \@div);
        @div = @{$_[0]};
        my @newfootprint = @{$_[1]};
        @footprint = map {$footprint[$_] > 0 || $newfootprint[$_]} (0.. $#div); 
        $rank = ($rank > $div[$i]) ? $div[$i] : $rank;
        unless ($rank) {
            $q = $i;
            last;
        }
    }
    return [$rank, $q];
}

sub _reduce {
=head1 reduce

Reduce div with respect to the q vertex

=cut
    my $q = shift;
    my @div = @{shift()};
    my $keep_intermediate_steps = shift;
    
    my @intermediate_steps;
    
    #Fireup the nodes that don't burn, while the graph doesn't burn
    while (my $burnt = _burn($q, \@div)) {
        #say "'", array2div(@div), "',";
        push @intermediate_steps, [@div] if ($keep_intermediate_steps); #we make a new reference [@div] 
    }
    
    return \@intermediate_steps if ($keep_intermediate_steps);
    return @div;
}

sub _reduce_footprint {
=head1 reduce

Reduce div with respect to the q vertex and keep track of the footprint left behind. Can also keep track of intermediate steps if desired.

=cut

    my $q = shift;
    my @div = @{shift()};
    my @footprint = map {$_ > 0} @div;
    my $keep_track_intermediate = shift;
    my @intermediate;
    
    #Fireup the nodes that don't burn, while the graph doesn't burn
    while (my $burnt = _burn($q, \@div)) {
        @footprint = map {$footprint[$_] > 0 || $div[$_]} (0.. $#div); 
        push @intermediate, array2div(@div) if ($keep_track_intermediate);
    }
    
    return \@div, \@footprint, \@intermediate;
}

sub _reduce_footprint_stop {
=head1 reduce

Reduce div with respect to the q vertex and keep track of the footprint left behind, stop when a chip falls to q

=cut

    my $q = shift;
    my @div = @{shift()}; 
    my @footprint = map {$_ > 0} @div;
    
    #Fireup the nodes that don't burn, while the graph doesn't burn
    while (my $burnt = _burn($q, \@div)) {
        last if ($div[$q]);
        @footprint = map {$footprint[$_] > 0 || $div[$_]} (0.. $#div); 
    }
    
    return \@div, \@footprint;
}



sub is_connected() {
=pod
is_connected(): Returns 1 if the Graph is connected, returns 0 if the Graph is not connected.
=cut    
    my $self = shift;
    
    $self->_context(); 
    my $burns = !_burn(0, $self->{empty_divisor}); 
    _clean_context();
    
    return $burns;
}

sub _burn {
=pod
_burn: PRIVATE Return false if the graph burns, otherwise modify in place the divisor and return a list of the nodes that burned 
=cut
    no warnings 'uninitialized';
    
    my $q = shift;
    my $div = shift();
    my (@burnt, @is_burnt, @burning, @fires);
    @burnt = ($q);
    @burning = ($q);
    @fires = @{$div};
    $fires[$q] = -1;
    $is_burnt[$q] = 1;

    #Spread fire
    while (scalar @burning > 0) {
        foreach my $i (@burning) {
            $fires[$_] = $fires[$_] + $lap[$i][$_] foreach (@{$support_lap[$i]}); 
        }
        @burning = grep {
            push @burnt, $_ if (!$is_burnt[$_] && $fires[$_] < 0 && !$is_burnt[$_]++); #this is how we avoid to waste several additions
        } map {@{$support_lap[$_]}} @burning; #here we replace the burning array for its neighbors
    }

    if (scalar @burnt == $n) {
        return 0;
    } else {
        $fires[$q] += 1 + $div->[$q];
        @{$div} = @fires; 
        return(\@burnt);
    }
}

sub neighborhood {
=pod
Return a list with the neighborhood of the vertex
=cut

    return @{$_[0]->{support}->[$_[1]]}
}

sub div_approximate {
=pod
Given a divisor D, and a number a, list all posibilities of approximating chips at vertices v_i to its neighbors, when i >= a. The code only works well when there is a single chip at a vertex, and that's the chip we wish to approximate. Multiple chips in different vertices is fine as long as we only have one chip per vertex. In practice this is the only situation we are currently interested in.
=cut

    my $self = shift;
    $self->_context();
    my $div = shift;
    my $a = shift;
    
    my @div_array = _div2array($div);
    
    #The whole setup
    $div =~ s/\d*v//;
    my @div = split / \d+v/, $div;
    @div = grep {$_ >= $a} @div;
    my @approx = (\@div);
    
    #We do the substitutions with a clever use of map ability to replace an element by various elements
    foreach my $i (0..$#div) {
        next if ($div[$i] < $a);
        
        $div_array[$div[$i]] = undef;
        @approx = map {my @copy = @{$_}; map {$copy[$i] = $_; [@copy]} $self->neighborhood($copy[$i])} @approx;
    }
    
    #Now we translate everything back to a divisor
    @approx = map {
        my @copy_div_array = @div_array;
        $copy_div_array[$_]++ foreach @{$_};
        array2div(@copy_div_array)
    } @approx;
    
    return @approx;
}

sub some_equivalent {
=pod
Given a list of vertices, reduce with respect to those vertices and keep track of the intermediate steps. If no list is given, reduce with respect to all vertices. 
=cut
    my $self = shift;
    $self->_context();
    my $div = shift;
    my @reduce = ($_[0]) ? @{$_[0]} : @vertices;
    my @equiv;
    
    #push foreach my $q (@reduce);
}

#=====================Energy Methods===================
sub calculate_Li {
    my $self = shift;
    $self->_context();
    my $q = shift; $q++; #minor method takes indices starting at 1
    
    my $mat = Math::MatrixReal->new_from_rows($self->{lap})->minor($q,$q);
    my @test;
    
    foreach my $i (0..$n-2) {
        foreach my $j (0..$i) {
            $test[$i][$j] = round($mat->minor($i+1,$j+1)->det());
            $test[$j][$i] = $test[$i][$j];
        }
    }
    print _matrix2string(@test);

    
    my $det = round($mat->det()); 
    $mat = $mat / $det;
    $mat = $mat->inverse();
    my @Li;
    
    foreach my $i (1..$self->{v}-1) {
        $Li[$i-1][$_-1]  = round($mat->element($i,$_)) foreach (1..$self->{v}-1);
        
    }
    
    $self->{det} = $det;
    $self->{Li} = \@Li;
    $self->{ginverses}{$q} = \@Li;
}

#============================These are presentation functions=========================


sub _div2array {
    my @div = @empty_divisor; 
    
    foreach (split(' ', $_[0])) {
        @_ = split 'v';
        $div[$_[1]] += $_[0];
    } 
    
    return @div;
}

sub array2div {
    my @div; 
    for (@vertices) {
        push @div, $_[$_] . 'v' . $_ if ($_[$_]);
    }
    
    return join(' ', @div);
}


sub _convert_div {
        my $ite = shift;
        my @div = @empty_divisor;
        my $j = 0;
        for my $i (@$ite) {
            if ($i) {
                $div[$j]++
            } else {
                $j++
            }
        }
        return @div;
}

sub _matrix2string {
    return join "\n", map {join "\t", @{$_}} @_; 
}

sub _string2matrix {
    return [map { [split "\t"] } split "\n", $_[0]];
}

sub _links2string {
    return join "\n", map {"$_->[1],$_->[2] "} @_;
}

sub _links2matrix {
    my @mat;
    my ($n, $i, $j);
    $n = 0;
    
    #convert to chain, this is a clever contraption
    map { ($i, $j) = split ',', $_; 
        $mat[$i][$j]++;
        $mat[$j][$i]++;
        $n = $j if ($n < $j);
        $n = $i if ($n < $i);} 
        split ' ' foreach (@_);
    
    #pad with zeros
    foreach my $ref (@mat) {
        @$ref = map {($ref->[$_]) ? $ref->[$_] : 0} (0..$n);
    }
    
    return \@mat;
}


#=======================================================================================
#                               Private Methods
#=======================================================================================
=pod
Iterators that generate all the possible divisors on a graph. 
I believe we wanted to do something different for the energy methods, to use a different iterator.
=cut

sub _all_divisors_iterator1 {
=pod
    We receive $c, the number of stones to use
    This is an optimized version for the purpose of skipping non-reduced divisors and finding all the rank1 divisors of a given graph. For the general purpose iterators see the next ones.
=cut

    #This is based in listing all possible ways of adding $c using $v numbers
    my $c = shift();
    my $s = $n + $c - 1;
    my $ite = [];
    push(@$ite, $_) for (1..$c);

    #Define the state and print the first one
    push(@$ite, (0) x ($s-$c));
    my $last = $c-1;
    my $previous_failure = 0;
    my $div_last = 0;
    my $is_finished;
    my $burns = 1;
    
    my @div = _convert_div($ite);
    $div[0] += 1;
    
    my @divisors;
    my $reduced_divisors;
    
    #begin timing 
    my $t0 = [gettimeofday];
    
    while (1) {
        if ($burns and !_burn(0, \@div)) {
            $reduced_divisors++;
            
            my $rank = _rank(\@div);
                if ($rank->[0] > 0) {
                    my $div_string = array2div(@div);
                    push @divisors, [$div_string, $rank->[0]];
                }
        } else {$burns = 1}
    
        last if ($is_finished);
        #Check if there is space to move $last to the next slot
        if ($last < ($s-1) ) {
            $ite->[$last] = 0;
            $last++;
            $ite->[$last] = $last + 1;
        } else {
            #Check how many elements we have at the end, begin at the end and break on the first zero we encounter
            my $sum = 0;
            for (1..$c) {
            last if !$ite->[$s - $_];
            $sum++;
            }

            #travel backwards from $last to find the previous 1 
            my $next;
            for (0..($last-$sum)) {
            $next = $_ if $ite->[$_];
            }
            #bump the last one 
            $ite->[$next] = 0;
            $ite->[$next+1] = $next+2;
            
            #and bring back what has accumulated at the end
            ##print "\n($last, $next)!\n";
            #It is two loops because otherwise we got a bug
            $ite->[$last + 1 - $_] = 0 for (1..$sum);
            $ite->[$next + 1 + $_] = $next + 2 + $_ for (1..$sum);

            $last = $next + 1 + $sum;
        }
        #convert the ite to div
        @div =  _convert_div($ite); 
        $is_finished = $div[-1] == $c;

        #convert the "last" place to a div place. Check if we fullfill the degree condition. If we failed in the previous one, check we don't fail now.
        $div[0] += 1;
        $div_last = $last - $c + 1;
        next if ($div[$previous_failure] >= $degrees[$previous_failure] and $previous_failure != 0 and !($burns = 0));
        if ($div[$div_last] >= $degrees[$div_last] and $div_last != 0 and !($burns = 0)) {
            $previous_failure = $div_last;
            next;
        }
        
    }
    return {rank1 => \@divisors, time_interval => tv_interval($t0), reduced_divisors => $reduced_divisors};
}


sub _all_divisors_iterator2 {
=pod
    We receive $c, the number of stones to use
    And $f, the function to call with the divisor
=cut

    #This is based in listing all possible ways of adding $c using $v numbers
    my $c = shift;
    my $test = shift;
    my @arguments = @_; 
    my $s = $n + $c - 1;
    my $ite = [];
    push(@$ite, $_) for (1..$c);

    #Define the state and print the first one
    push(@$ite, (0) x ($s-$c));
    my $last = $c-1;
    my $previous_failure = 0;
    my $div_last = 0;
    my $is_finished;
    my @div = _convert_div($ite);
    my $burns;
    
    my @divisors;
    my $reduced_divisors;
    
    #begin timing 
    my $t0 = [gettimeofday];
    
    while (1) {
        push @divisors, $_ if ($_ = &$test(\@div, @arguments));
        last if ($is_finished);
        
        #Check if there is space to move $last to the next slot
        if ($last < ($s-1) ) {
            $ite->[$last] = 0;
            $last++;
            $ite->[$last] = $last + 1;
        } else {
            #Check how many elements we have at the end, begin at the end and break on the first zero we encounter
            my $sum = 0;
            for (1..$c) {
            last if !$ite->[$s - $_];
            $sum++;
            }

            #travel backwards from $last to find the previous 1 
            my $next;
            for (0..($last-$sum)) {
            $next = $_ if $ite->[$_];
            }
            #bump the last one 
            $ite->[$next] = 0;
            $ite->[$next+1] = $next+2;
            
            #and bring back what has accumulated at the end
            ##print "\n($last, $next)!\n";
            #It is two loops because otherwise we got a bug
            $ite->[$last + 1 - $_] = 0 for (1..$sum);
            $ite->[$next + 1 + $_] = $next + 2 + $_ for (1..$sum);

            $last = $next + 1 + $sum;
        }
        #convert the ite to div
        @div =  _convert_div($ite); 
        $is_finished = $div[-1] == $c;
    }
    return @divisors;
}

=pod
To ease the names of variables in the code and do some code speedups, the routines here work by having a context.
Essentially the routines in the API call the _context() private method and then proceed to pretend that the global 
variables are safe and happy. Therefore if you use a Private Method, be sure that the context is properly set up.
=cut

sub _context() {
    my $self = shift;
    $n = $self->{n};
    @mat = @{$self->{matrix}};
    @lap = @{$self->{lap}};
    @support = @{$self->{support}};
    @support_lap = @{$self->{support_lap}};
    @vertices = @{$self->{vertices}};
    @vertices_sorted = @{$self->{vertices_sorted}};
    @empty_divisor = @{$self->{empty_divisor}};
    @degrees = @{$self->{degrees}};
    $c = $self->{c};
}

sub _clean_context() {
    $n = undef;
    @mat = undef;
    @lap = undef;
    @support = undef;
    @support_lap = undef;
    @vertices = undef;
    @vertices_sorted = undef;
    @empty_divisor = undef;
    @degrees = undef;
    $c = undef;
}

sub use_database {
=head1 use_database
There are several calculations that take a while, so we save them in a database
=cut
    my $db_file = $_[1];
    $table_name = $_[2];
    
    $dbh = DBI->connect("dbi:SQLite:dbname=$db_file")
        or die $DBI::errstr;
    
    #We have different tables for different experiments
    #Table to save all the rank1 divisors
    $dbh->do("CREATE TABLE IF NOT EXISTS $table_name\_rank1 (graph TEXT, divisors TEXT, rank1 INT, rank2 INT, total INT, reduced INT, g INT, v INT, signature TEXT, det INT, factorization TEXT, time INT)" );
    $dbh->do("CREATE INDEX IF NOT EXISTS $table_name\_total ON $table_name\_rank1 (total)" );
    
    #Table to save the count of reduced divisors
}

sub use_table {
    $table_name = $_[1];
    
    #We have different tables for different experiments
    #Table to save all the rank1 divisors
    $dbh->do("CREATE TABLE IF NOT EXISTS $table_name\_rank1 (graph TEXT, divisors TEXT, rank1 INT, rank2 INT, total INT, reduced INT, g INT, v INT, signature TEXT, det INT, factorization TEXT, time INT)" );
    $dbh->do("CREATE INDEX IF NOT EXISTS $table_name\_total ON $table_name\_rank1 (total)" );
}

sub commit_database {
=head1 use_database
Given a table name and a value string commit the values to the database
=cut
    my $suffix = $_[1];
    my $values = $_[2];

    $dbh->do("INSERT INTO $table_name\_$suffix VALUES $values");
}

sub BUILD {
    #We create a default database connection unless there is a pre-existing database connection already
    ChipFiring->use_database('graphs.sql', 'graphs') unless ($dbh);
}

1;
