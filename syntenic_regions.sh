#!/bin/sh -e

# Note on input: one-to-one optimal alignment with last-split in both
# directions.  In principle, running maf-swap before or after last-split should
# not matter, but in corner cases where two alternate alignments are equally
# good, the following orders ensures that the chosen scaffold is the one sorted
# first in the query genome file (in our case it means it is longer.
# maf-swap | last-split | maf-swap | last-split

# Add the target sequence name as first field for sorting.
# The regular expression assumes scaffold names do not contain dots.
perl -pe '
  next if /\#/ ;
  ($target) = /ID=(.+?)[\.]/ if /\tmatch\t/ ;
  ($target) = /Parent=(.+?)[\.:]/ if /\tmatch_part\t/ ;
  $_ = join "\t", $target, $_
' |
# Sort on first field, then seqname, then start position.
# The reverse flag on the first key ensures the header stays at the top.
# To also sort on strand add -k8,8 after the first sort key/
sort -k1,1Vr -k2,2V -k5,5n |
# add a Parent attribute to the consecutive lines that are
# on the same seqname, same strand and within 500,000 nt of each other.
perl -ape '
  use constant {target => 0, seqname => 1, type => 3, start => 4, end => 5};
  next if /\#/ ;
  s/.+?\t// ; # Remove sort key
  next unless $F[type] eq "match" ;
  chomp ;
  $parent = "$F[target]→$F[seqname]:$F[start]"
    unless ($P[target]      eq $F[target]  &&
            $P[seqname]     eq $F[seqname] &&
            $P[end] + 500000 > $F[start]) ;
  @P = @F ;
  $_ = "$_;Parent=$parent\n"
' |
# Add a dummy line so that the last syntenic region is also computed.
sed '$askipme\t.\tmatch\t.\t.\t.\t.\t.\t.\tParent=skipme' |
# Output a "syntenic_region" line grouping each line that have the same Parent.
# Swap commented line for stranded output (but also see the sort command above).
perl -apE '
  use constant {seqname => 0, type => 2, start => 3, end => 4, strand => 6};
  next if /\#/ ;
  next unless $F[type] eq "match" ;
  ($curparent) = /Parent=(.+)$/ ;
  if (defined ($parent) && $parent ne $curparent) {
    $name = $parent =~ s/→.*$//r ;
    #say join("\t", $P[seqname], ".", "syntenic_region", $P[start], $P[end], ".", $P[strand], ".", "ID=$parent;Name=$name") ;
    say join("\t", $P[seqname], ".", "syntenic_region", $P[start], $P[end], ".", ".", ".", "ID=$parent;Name=$name") ;
    @P = @F ;
    last if $curparent eq "skipme" ;
    $parent = $curparent ;
    next} ;
  @P = @F if not @P ;
  $P[start] = $F[start] if $F[start] < $P[start] ;
  $P[end] = $F[end] if $F[end] > $P[end] ;
  $parent = $curparent'
