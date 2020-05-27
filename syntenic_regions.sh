#!/bin/sh -e

# One-to-one optimal alignment with last-split in both orders.
# Make sure maf-swap is not the last command of the pipe.
last-split -m1 $1 | maf-swap | last-split -m1 | maf-swap |
# Convert to GFF3, join blocks within 200,000 nt.
~/maf-convert gff -j 200000 |
# Add the target sequence name as first field.
# The regular expression assumes scaffold names do not contain dots.
perl -pe '
  next if /\#/ ;
  ($target) = /ID=(.+?)[\.]/ if /\tmatch\t/ ;
  ($target) = /Parent=(.+?)[\.:]/ if /\tmatch_part\t/ ;
  $_ = join "\t", $target, $_
' |
# Sort on strand, then first field, then seqname, then start position.
# The reverse flag on the first key ensures the header stays at the top.
sort -k8,8 -k1,1Vr -k2,2V -k5,5n |
# add a Parent attribute to the consecutive lines that are
# on the same seqname, same strand and within 500,000 nt of each other.
perl -ape '
  use constant {target => 0, seqname => 1, type => 3, start => 4, end => 5};
  next if /\#/ ;
  s/.+?\t// ; # Remove sort key
  next unless $F[type] eq "match" ;
  chomp ;
  $parent = "$F[target]_$F[start]"
    unless ($P[target]      eq $F[target]  &&
            $P[seqname]     eq $F[seqname] &&
            $P[end] + 500000 > $F[start]) ;
  @P = @F ;
  $_ = "$_;Parent=$parent\n"
' |
# Output a "syntenic_region" line grouping each line that have the same Parent.
perl -apE '
  use constant {seqname => 0, type => 2, start => 3, end => 4, strand => 6};
  next if /\#/ ;
  next unless $F[type] eq "match" ;
  ($curparent) = /Parent=(.+)$/ ;
  if (defined ($parent) && $parent ne $curparent) {
    $name = $parent =~ s/_[^_]*$//r ;
    say join("\t", $P[seqname], ".", "syntenic_region", $P[start], $P[end], ".", $P[strand], ".", "ID=$parent;Name=$name") ;
    @P = @F ;
    $parent = $curparent ;
    next} ;
  @P = @F if not @P ;
  $P[start] = $F[start] if $F[start] < $P[start] ;
  $P[end] = $F[end] if $F[end] > $P[end] ;
  $parent = $curparent'
