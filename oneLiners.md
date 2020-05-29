Total length of short tandem repeats
```
$ awk '{print $3 - $2}' I69-4.tantan.f4.txt | paste -sd+ - | bc
1497819
```

ACGT composition
```
module load emboss; $ compseq -word 1 -filter I69-4.fa
[…]
A	18934744		0.2940680	0.2500000	1.1762720
C	13210172		0.2051619	0.2500000	0.8206478
G	13221178		0.2053329	0.2500000	0.8213315
T	19004406		0.2951499	0.2500000	1.1805995
```
