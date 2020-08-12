#!/usr/bin/awk -f

BEGIN {
  OFS="\t"
  print "##gff-version 3"
}

{
  print $1, "tantan", "tandem_repeat", $2 + 1, $3, int($4*$5), ".", ".", "Name="$6
}
