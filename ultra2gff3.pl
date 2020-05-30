#!/usr/bin/perl -w

use feature "say" ;
use constant {seqid  => 0, source => 1, type  => 2,
              start  => 3, end    => 4,
              score  => 5, strand => 6, phase => 7, attrs => 8} ;

my @line                        ;
$line[source] = "ultra2gff3.pl" ;
$line[type]   = "tandem_repeat" ;
$line[phase]  = "."             ;
$line[strand] = "+"             ;

my $headerIsOver = 0 ;

say "##gff-version 3" ;

while (<>) {
  /PassID/                  and $headerIsOver = 1                         ;
  next unless $headerIsOver                                               ;
  /SequenceName.+?"(\w+)"/  and $line[seqid]  = $1                        ;
  /Start.+?(\d+)/           and $line[start]  = 1 + $1                    ;
  /Length.+?(\d+)/          and $line[end]    = $line[start] + $1         ;
  /Period.+?(\d+)/          and $line[attrs]  = "Period=$1;" ; # Reset attrs
  /Score.+?([\d\.]+)/       and $line[score]  = $1                        ;
  /Log2 Pval.+?([\d\.+-]+)/ and $line[attrs] .= "Log2Pval=$1;"            ;
  /Substitutions.+?(\d+)/   and $line[attrs] .= "Substitutions=$1;"       ;
  /Insertions.+?(\d+)/      and $line[attrs] .= "Insertions=$1;"          ;
  /Deletions.+?(\d+)/       and $line[attrs] .= "Deletions=$1;"           ;
  /Consensus.+?"(\w+)"/     and $line[attrs] .= "Consensus=$1;Name=$1;"   ;
#  /Sequence.+?"(\w+)"/      and $line[attrs] .= "Sequence=$1;"            ;
  /}/                       and say join "\t", @line
}
