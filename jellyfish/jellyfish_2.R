#!/bin/Rscript

library(vioplot)

genomescope_data <- read.delim('GenomeScope_Jellyfish_2_JustOki.csv',header=T, stringsAsFactors = F)

# The full data set contains multiple sequencing runs
# as well as non-converged models, as well as models with
# aberrant results. 
#
# The subset from 57-63, which is based on the MiSeq
# polishing reads, varies across a few values for K 
# but excludes the worst K values (13 or >41).
genomescope_of_interest <- genomescope_data[c(57:63),]

gs <- as.data.frame(cbind(genomescope_of_interest$genome_repeat_min, genomescope_of_interest$genome_repeat_max, genomescope_of_interest$genome_unique_min, genomescope_of_interest$genome_unique_max, genomescope_of_interest$genome_haploid_min, genomescope_of_interest$genome_haploid_max))

gs <- apply(gs, 2, function(x) x/1000000)

rownames(gs) <- genomescope_of_interest$k
colnames(gs) <- c('Repetitive genome (min)','Repetitive genome (max)', 'Unique genome (min)', 'Unique genome (max)', 'Haploid genome (min)', 'Haploid genome (max)')


pal.3 <- c('#00d87f', '#ce4dd0', '#00b5c8')
pal.6 <- c('#a6cee3','#1f78b4', '#b2df8a', '#33a02c','#fb9a99', '#e31a1c')

# base R boxplot
p1 <- function() {
  par(mar=c(5,12,2,2))
  boxplot(gs, horizontal=TRUE, xlab='Genomescope estimated size (Mb)', col=pal.6, pch=21, bg='white', yaxt='n')
  axis(2, at=1:6, labels=colnames(gs), las=2)
  abline(v=c(20,30,40,50), lty=3, col='#d3d3d3')
}

#same thing but violin plot
p2 <- function() {
  par(mar=c(5,12,2,2))
  vioplot(gs[,1], gs[,2], gs[,3], gs[,4], gs[,5], gs[,6], col=pal.6, horizontal=TRUE, yaxt='n', xlab='Genomescope estimated size (Mb)')
  axis(1, at = c(20, 30, 40, 50), labels=c(20,30,40,50) )
  axis(2, at = c(1:6), labels=colnames(gs), las=2 )
  abline(v=c(20,30,40,50), lty=3, col='#d3d3d3')
}

# The min and max are so close to one another that
# it might make sense to just combine them into 
# a single distribution instead.
gs2 <- t(as.data.frame(rbind(c(gs[,1], gs[,2]), c(gs[,3], gs[,4]), c(gs[,5], (gs[,6])))))
colnames(gs2) <- c('Repetitive genome', 'Unique genome', 'Haploid genome')

p3 <- function() {
  boxplot(gs2, horizontal=TRUE, xlab='Genomescope estimated size (Mb)', col=pal.3, pch=21, bg='white', yaxt='n', main='')
  axis(2, at=1:3, labels=colnames(gs2), las=2)
  abline(v=c(20,30,40,50), lty=3, col='#d3d3d3')
}

p4 <- function() {
  par(mar=c(5,12,2,2))
  vioplot(gs2[,1], gs2[,2], gs2[,3], col=pal.3, horizontal=TRUE, yaxt='n', xlab='Genomescope estimated size (Mb)')
  axis(1, at = c(20, 30, 40, 50), labels=c(20,30,40,50) )
  axis(2, at = c(1:3), labels=colnames(gs2), las=2 )
  abline(v=c(20,30,40,50), lty=3, col='#d3d3d3')
}


pdf('jellyfish_2.pdf')
p1()
p2()
p3()
p4()
dev.off()