#!/bin/Rscript

library(circlize)

genome_size = 9225
mito <- read.delim('HiCscaffold40.gff', header=F, stringsAsFactors = F)
# Remove empty column at the end of the data frame
mito <- mito[,-10]
# rename columns according to gff3 specifications
colnames(mito) <- c('sequence','source','feature','start','end','score','strand','phase','attributes')

# Complicated string splitting to get just gene names
mito$name <- unname(gsub(x=sapply(mito$attributes, function(x) unlist(strsplit(x, '='))[2] ), pattern='\\_', replacement=' '))

# Gene midpoints for text labels
mito$midpoint <- (mito$start+mito$end)/2

##############################################
#
# Plotting parameters. Global options up front.
#
##############################################
pal.10 <- c('#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69','#fccde5','#bc80bd')

mito$colour <- pal.10[10]
mito$colour[grep('RNA', mito$feature)] <- pal.10[4]
mito$colour[grep('gene', mito$feature)] <- pal.10[5]
mito$colour[grep('rep_origin', mito$feature)] <- pal.10[2]

# Easier to split the data frame into positive and 
# negative strands.
mito_plus <- subset(mito, strand == '+')
mito_minus <- subset(mito, strand == '-')

pdf('Mitochondrion.pdf')
# Initialization for track sizes
circos.par(track.height=0.10)
circos.initialize(c('1'), xlim=c(0,genome_size))
circos.track(ylim=c(0,1)) # present or absent basically
circos.axis(h=1.60, major.tick=1000, labels=FALSE)
circos.axis(h=1.60, major.at=c(0, 2000, 4000, 6000, 8000, 9225), labels=c('0kb', '2kb', '4kb', '6kb', '8kb', '9.225kb'), major.tick=F)

# Sense strand first.
circos.rect(xleft=mito_plus$start, xright=mito_plus$end, ybottom=rep(0, nrow(mito_plus)), ytop=rep(1, nrow(mito_plus)), col=mito_plus$colour, border='black')

circos.text(x=mito_plus$midpoint, y=rep(-0.75,nrow(mito_plus)), mito_plus$name, cex=0.75, niceFacing = T)

# Antisense strand on the inside to include more labels.
circos.par(track.margin=c(0.125,0.125))
circos.track(ylim=c(0,1))
circos.rect(xleft=mito_minus$start, xright=mito_minus$end, ybottom=rep(0, nrow(mito_minus)), ytop=rep(1, nrow(mito_minus)), col=mito_minus$colour, border='black')
circos.text(x=mito_minus$midpoint, y=rep(-0.75,nrow(mito_minus)), mito_minus$name, cex=0.75, niceFacing = T)

legend('topright', pch=15, col=c(pal.10[5],pal.10[4],pal.10[2]), legend=c('Gene','RNA','Origin of replication'))

dev.off()
# Make an antisense track
circos.clear()