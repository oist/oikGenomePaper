#!/bin/Rscript

library(splines)

lengths_genes <- read.delim('length_genes.tsv', header=T, stringsAsFactors = F)
lengths_exons <- read.delim('length_exons.tsv', header=T, stringsAsFactors = F)
lengths_introns <- read.delim('length_introns.tsv', header=T, stringsAsFactors = F)

lengths_genes[,1] <- as.numeric(lengths_genes[,1])
lengths_genes[,2] <- as.numeric(lengths_genes[,2])
lengths_exons[,1] <- as.numeric(lengths_exons[,1])
lengths_exons[,2] <- as.numeric(lengths_exons[,2])
lengths_introns[,1] <- as.numeric(lengths_introns[,1])
lengths_introns[,2] <- as.numeric(lengths_introns[,2])

lengths_genes$probability <- as.numeric(sapply(lengths_genes$distribution, function(x) unname(unlist(strsplit(unlist(strsplit(x, '='))[2], ','))[1]) ))
lengths_exons$probability <- as.numeric(sapply(lengths_exons$distribution, function(x) unname(unlist(strsplit(unlist(strsplit(x, '='))[2], ','))[1]) ))
lengths_introns$probability <- as.numeric(sapply(lengths_introns$distribution, function(x) unname(unlist(strsplit(unlist(strsplit(x, '='))[2], ','))[1]) ))

maxprob_gene <- lengths_genes[which.max(lengths_genes$probability),1]
maxprob_exon <- lengths_exons[which.max(lengths_exons$probability),1]
maxprob_intron <- lengths_introns[which.max(lengths_introns$probability),1]

pdf('transcriptome_properties.pdf')
par(mfrow=c(3,1))

plot(lengths_genes$length~lengths_genes$gene, type='p', col=adjustcolor('#fb9a99', 0.1), pch=16, ylab='Count', xlab='Length', main='Gene length distribution (bp)')
lines( smooth.spline(lengths_genes$length~lengths_genes$gene), col='#e41a1c')
abline(v=maxprob_gene, col='#595959', lty=3)
text(x=10000, y=12, labels=paste("Maximum prob.:", maxprob_gene, 'bp'))

plot(lengths_exons$length~lengths_exons$exon, type='p', col=adjustcolor('#a6cee3', 0.1), pch=16,  ylab='Count', xlab='Length', main='Exon length distribution (bp)')
lines( smooth.spline(lengths_exons$length~lengths_exons$exon), col='#1f78b4')
abline(v=maxprob_exon, col='#595959', lty=3)
text(x=3000, y=300, labels=paste("Maximum prob.:", maxprob_exon, 'bp'))

plot(lengths_introns$length~lengths_introns$intron, type='p', col=adjustcolor('#b2df8a', 0.1), pch=16,  ylab='Count', xlab='Length', main='Intron length distribution (bp)')
lines( smooth.spline(lengths_introns$length~lengths_introns$intron), col='#33a02c')
abline(v=maxprob_intron, col='#595959', lty=3)
text(x=10000, y=3000, labels=paste("Maximum prob.:", maxprob_intron, 'bp'))

plot( lengths_genes$length~log10(lengths_genes$gene), type='p', col=adjustcolor('#fb9a99', 0.1), pch=16, ylab='Count', xlab='log10(Length (bp))', main='log10(Gene length distribution)')
lines( smooth.spline( lengths_genes$length~log10(lengths_genes$gene)), col='#e41a1c')
abline(v=log10(maxprob_gene), col='#595959', lty=3)
text(x=3.5, y=15, labels=paste("Maximum prob.:", maxprob_gene, 'bp'))

plot(lengths_exons$length~log10(lengths_exons$exon), type='p', col=adjustcolor('#a6cee3', 0.1), pch=16,  ylab='Count', xlab='log10(Length (bp))', main='log10(Exon length distribution)')
lines( smooth.spline(lengths_exons$length~log10(lengths_exons$exon)), col='#1f78b4')
abline(v=log10(maxprob_exon), col='#595959', lty=3)
text(x=3, y=300, labels=paste("Maximum prob.:", maxprob_exon, 'bp'))

plot(lengths_introns$length~log10(lengths_introns$intron), type='p', col=adjustcolor('#b2df8a', 0.1), pch=16,  ylab='Count', xlab='log10(Length (bp))', main='log10(Intron length distribution)', ylim=c(-500,6335))
lines( smooth.spline(lengths_introns$length~log10(lengths_introns$intron)), col='#33a02c')
abline(v=log10(maxprob_intron), col='#595959', lty=3)
text(x=2.25, y=3000, labels=paste("Maximum prob.:", maxprob_intron, 'bp'))


dev.off()
