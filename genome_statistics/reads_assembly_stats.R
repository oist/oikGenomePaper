#!/bin/Rscript

library(data.table)

distribution_reads <- fread('read_length_distribution.gz',header=F, stringsAsFactors = F)
# the Contig_Coverage.tsv file has lengths and coverage values for every contig
longread_coverage <- fread('Contig_Coverage.tsv', header=T, stringsAsFactors = F)

distribution_reads <- distribution_reads[order(distribution_reads$V1, decreasing=T),]



# Reorder by contig length
longread_coverage <- longread_coverage[order(as.numeric(longread_coverage$Length), decreasing=T),]

# Transform contig lengths into megabases
longread_coverage$Length_Mb <- longread_coverage$Length/1000000
longread_coverage$Length_log10 <- log10(longread_coverage$Length)

# Calculating N50 from a numeric vector containing a distribution of lengths.
# This works for assemblies and read length distributions.
calc_N50 <- function(numeric_vector_of_lengths){
  tmp_df <- as.data.frame(matrix(nrow=length(numeric_vector_of_lengths), ncol=2))
  tmp_df$V1 <- numeric_vector_of_lengths
  tmp_df <- tmp_df[order(tmp_df$V1, decreasing=T),]
  
  tmp_df$V2 <- cumsum(as.numeric(tmp_df$V1))
  half_total_size <- 0.5 * sum(tmp_df$V1)
  N50 <- head(tmp_df[which(tmp_df$V2 > half_total_size),1], n=1)
  return(N50)
}

n50_assembly <- calc_N50(longread_coverage$Length)
n50_reads <- calc_N50(distribution_reads$V1)

############################################################
#
# Plotting commands, not just 
#
#
############################################################

pal.coldyellow_hotred <- colorRampPalette(c('#ffffcc','#ffeda0','#fed976','#feb24c','#fd8d3c','#fc4e2a','#e31a1c','#bd0026','#800026'))(10) 

pdf('reads_and_assembly_distribution.pdf')
# Read length distribution histogram
hist(distribution_reads$V1/1000000, xlab='Length (Mb)', col='lightblue', main='Distribution of read lengths')
# Pseudo-N50. Calculated in R. Sort by length, calculate cumulative sum, then
# N50 = value of X which is greater than or equal to 1/2 of cumulative sum.
abline(v=n50_reads/1000000, lty=3, col='#595959')
text(x=0.03, y=4000000, labels = paste('Pseudo-N50:', n50_reads) )

hist(longread_coverage$Length_Mb, xlab='Length (Mb)', col='pink', main='Distribution of contig lengths' )
abline(v=n50_assembly/1000000, lty=3, col='#595959')
text(x=13.5, y=30, labels=paste('N50:', round(n50_assembly/1000000, 3), 'Mb') )

# Untransformed plot
barplot(longread_coverage$Length_Mb, col=pal.coldyellow_hotred[ceiling(longread_coverage$`Sequana Depth of Coverage`/max(longread_coverage$`Sequana Depth of Coverage`)*10)], xlab='Contig (n=54)', ylab='Length (Mbp)', las=2, main='Contigs with coverage indicated')
abline(h=c(0,5,10,15), lty=3, col='#595959')
abline(v=2, lty=3, col='#595959')
text(x=12, y=16, labels = paste('N50:', round(n50_assembly/1000000, 3), 'Mb' ))
legend('topright', legend=c('~240X', '~7500X', '>17000X'), pch=15, col=c(pal.coldyellow_hotred[1], pal.coldyellow_hotred[5], pal.coldyellow_hotred[10]))

# Log-transformed plot
barplot(longread_coverage$Length_log10, col=pal.coldyellow_hotred[ceiling(longread_coverage$`Sequana Depth of Coverage`/max(longread_coverage$`Sequana Depth of Coverage`)*10)], xlab='Contig (n=54)', ylab='log10(Contig Length)', las=2, main='Contigs with log-transformed length and coverage indicated')
abline(h=c(0:7), lty=3, col='#595959')
legend('topright', legend=c('~240X', '~7500X', '~17000X'), pch=15, col=c(pal.coldyellow_hotred[1], pal.coldyellow_hotred[5], pal.coldyellow_hotred[10]))
abline(v=2, lty=3, col='#595959')
text(x=15, y=6, labels=paste('N50:', round(n50_assembly/1000000, digits=3), 'Mb') )


dev.off()
