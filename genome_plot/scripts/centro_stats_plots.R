#!/bin/Rscript

library(data.table)
library(parallel)

# NOTE! All of the below comes from the genome_plot.R script, until the next NOTE!
# It's copied and pasted here to get the inputs into a new format.

# Gather input options
window_file <- "../inputs/I69-5.windowed_stats.txt"
gene_file <- "../inputs/I69-5.genes.gff3"
gap_file <- "../inputs/I69-5.gaps.gff3"
output_file <- "Rplots.pdf"
window_size <- 50000
plot_replace_zeroes <- FALSE
make_cluster <- TRUE
num_threads <- 8

####################################################################
#
# Actual reading and writing of files for all of the work.
#
####################################################################

scaffold_sizes <- data.frame(scaffold=c('chr1', 'chr2', 'PAR', 'XSR', 'YSR'), chr_name=c('chr 1', 'chr 2', 'PAR', 'XSR', 'YSR'), starts=rep(1, 5), ends=c(14533022, 16158756, 17092476, 12959145, 2916375))

# Genome assembly parameters (GC and coverage over a window)
genome_stats <- fread(window_file, header=T, stringsAsFactors = F)
genome_stats$Window_Start <- as.numeric(sapply(genome_stats$Window, function(x) unlist(strsplit(x, '-'))[1]))
genome_stats$Window_End <- as.numeric(sapply(genome_stats$Window, function(x) unlist(strsplit(x, '-'))[2]))
genome_stats$Window_Mid <- (genome_stats$Window_Start + genome_stats$Window_End)/2
# Subset the assembly stats into the scaffolds of interest
genome_stats <- subset(genome_stats, Contig %in% scaffold_sizes$scaffold)
genome_stats$chr_name <- scaffold_sizes$chr_name[match(genome_stats$Contig, scaffold_sizes$scaffold)]
# Make sure the demoninator is correct depending on the window size of the files.
genome_stats$PctMasked <- round(genome_stats$NumMasked/window_size*100, digits=2)

# Gene data from gff file
genes <- fread(gene_file, stringsAsFactors=F, header=F, sep='\t', skip='#')
colnames(genes) <- c('sequence','source','feature','start','end','score','strand','phase','attributes')
# Subsetting the data to include only the chromosomal scaffolds
genes <- subset(genes, sequence %in% scaffold_sizes$scaffold)
genes <- subset(genes, feature == "gene")

# Rename the scaffolds according to chromosome names
genes$midpoint <- (genes$start+genes$end)/2
genes$chr_name <- scaffold_sizes$chr_name[match(genes$sequence, scaffold_sizes$scaffold)]

# Gaps
gaps <- fread(gap_file, stringsAsFactors=F, header=F, sep='\t', skip='#')
colnames(gaps) <- c('sequence','source','feature','start','end','score','strand','phase','attributes')
gaps$chr_name <- scaffold_sizes$chr_name[match(gaps$sequence, scaffold_sizes$scaffold)]

# Splitting into lists of chromosomes to ease plotting
split_genes <- split(x = genes, f = genes$chr_name)
split_gaps <- split(x = gaps, f = gaps$chr_name)
split_stats <- split(x = genome_stats, f = genome_stats$chr_name)

match_stat_coords_to_genes <- function(chr_name){
	chr_stats <- split_stats[chr_name][[1]]
	chr_genes <- split_genes[chr_name][[1]]
	gene_counts <- as.numeric(apply(chr_stats, 1, function(x) nrow(subset(chr_genes, midpoint > as.numeric(x[7]) & midpoint < as.numeric(x[8])))))
	chr_stats$gene_counts <- gene_counts
	return(chr_stats)
}

# Get gene densities for every scaffold, for every 5kb window
if(make_cluster){
	cluster <- makeForkCluster(num_threads)
	modified_split_stats <- parLapply(cl=cluster, X=names(split_stats), fun=match_stat_coords_to_genes)
	names(modified_split_stats) <- names(split_stats)
	split_stats = modified_split_stats
	rm(modified_split_stats)
	stopCluster(cluster)
	} else {
	modified_split_stats <- lapply(names(split_stats), function(x) match_stat_coords_to_genes(x))
	split_stats <- modified_split_stats
	rm(modified_split_stats)
}

# Calculate overall spline minimums and maximums.
# These are needed if smoothed_spline=TRUE and scaled_per_chromosome=FALSE.
overall_gc_min <- min(unlist(lapply(split_stats, function(x) smooth.spline(x$GC~x$Window_Mid)$y)))
overall_gc_max <- max(unlist(lapply(split_stats, function(x) smooth.spline(x$GC~x$Window_Mid)$y)))
overall_cov_min <- min(unlist(lapply(split_stats, function(x) smooth.spline(x$Depth~x$Window_Mid)$y)))
overall_cov_max <- max(unlist(lapply(split_stats, function(x) smooth.spline(x$Depth~x$Window_Mid)$y)))
overall_rep_min <- min(unlist(lapply(split_stats, function(x) smooth.spline(x$PctMasked~x$Window_Mid)$y)))
overall_rep_max <- max(unlist(lapply(split_stats, function(x) smooth.spline(x$PctMasked~x$Window_Mid)$y)))
overall_gene_min <- min(unlist(lapply(split_stats, function(x) smooth.spline(x$gene_counts~x$Window_Mid)$y)))
overall_gene_max <- max(unlist(lapply(split_stats, function(x) smooth.spline(x$gene_counts~x$Window_Mid)$y)))


# Special options for a global minimum of 0, ignoring whatever the spline was doing.
if(plot_replace_zeroes){
	if(overall_gc_min < 0){
		overall_gc_min <- 0
	}
	if(overall_cov_min < 0){
		overall_cov_min <- 0
	}
	if(overall_rep_min < 0){
		overall_cov_min <- 0
	}
	if(overall_gene_min < 0){
		overall_gene_min <- 0
	}
}
# A function to turn data into a smoothed spline, with the option to replace negative values with 0.
spline_coords <- function(x, y, makeNonZero=FALSE){
	coords <- smooth.spline(y~x)$y
	if(makeNonZero==TRUE){
		coords[coords < 0] <- 0
	}
	return(coords)
}


# NOTE! This is the end of the genome_plot.R script.
####################################################################
#
# The actual statistical stuff goes here.
#
####################################################################
scaffold_sizes$centromere_start <- c(5191657,5707009,6029625,NA,NA)
scaffold_sizes$centromere_end <- c(5192156,5707508,6030124,NA,NA)

chr_1_short <- subset(split_stats[['chr 1']], Window_End < scaffold_sizes[1,]$centromere_start )
chr_2_short <- subset(split_stats[['chr 2']], Window_End < scaffold_sizes[2,]$centromere_start )
chr_PAR_short <- subset(split_stats[['PAR']], Window_End < scaffold_sizes[3,]$centromere_start )
chr_1_long <- subset(split_stats[['chr 1']], Window_End > scaffold_sizes[1,]$centromere_end )
chr_2_long <- subset(split_stats[['chr 2']], Window_End > scaffold_sizes[2,]$centromere_end )
chr_PAR_long <- subset(split_stats[['PAR']], Window_End > scaffold_sizes[3,]$centromere_end )


# Doing some statistics...
# First, merge and reformat the data into a nice combined data frame
chr_1_short$chr_arm <- 'short'
chr_2_short$chr_arm <- 'short'
chr_PAR_short$chr_arm <- 'short'
chr_1_long$chr_arm <- 'long'
chr_2_long$chr_arm <- 'long'
chr_PAR_long$chr_arm <- 'long'
chr_1_j <- as.data.frame(rbind(chr_1_short, chr_1_long))
chr_2_j <- as.data.frame(rbind(chr_2_short, chr_2_long))
chr_PAR_j <- as.data.frame(rbind(chr_PAR_short, chr_PAR_long))
chr_list <- list('chr 1'=chr_1_j, 'chr 2'=chr_2_j, 'PAR'=chr_PAR_j)

test_significant_difference <- function(chr_name, test_function){
	chr_data <- chr_list[[chr_name]]
	if(as.character(substitute(test_function)) == 'ks.test'){
		short_arm <- subset(chr_data, chr_arm == 'short')
		long_arm  <- subset(chr_data, chr_arm == 'long')
		test_gc <- ks.test(short_arm$GC, long_arm$GC)
		test_depth <- ks.test(short_arm$Depth, long_arm$Depth)
		test_rep <- ks.test(short_arm$PctMasked, long_arm$PctMasked)
		test_genes <- ks.test(short_arm$gene_counts, long_arm$gene_counts)
		test_rs <- ks.test(short_arm$RESites, long_arm$RESites)

		tests <- list(test_gc, test_depth, test_rep, test_genes, test_rs)

		test_result_df <- as.data.frame(t(data.frame(c(test_gc$data.name, test_gc$method, test_gc$statistic, test_gc$p.value, test_gc$alternative), c(test_depth$data.name, test_depth$method, test_depth$statistic, test_depth$p.value, test_depth$alternative), c(test_rep$data.name, test_rep$method, test_rep$statistic, test_rep$p.value, test_rep$alternative), c(test_genes$data.name, test_genes$method, test_genes$statistic, test_genes$p.value,test_genes$alternative), c(test_rs$data.name, test_rs$method, test_rs$statistic, test_rs$p.value, test_rs$alternative))))

		rownames(test_result_df) <- NULL
		test_result_df <- cbind(chr_name, test_result_df)
		colnames(test_result_df) <- c('chr_name', 'formula', 'test', 'statistic', 'p_value', 'alternative_hypothesis')

	} else {
		test_gc <- test_function(chr_data$GC~chr_data$chr_arm)
		test_depth <- test_function(chr_data$Depth~chr_data$chr_arm)
		test_rep <- test_function(chr_data$PctMasked~chr_data$chr_arm)
		test_genes <- test_function(chr_data$gene_counts~chr_data$chr_arm)
		test_rs <- test_function(chr_data$RESites~chr_data$chr_arm)

		tests <- list(test_gc, test_depth, test_rep, test_genes, test_rs)

		test_result_df <- as.data.frame(t(data.frame(c(test_gc$data.name, test_gc$method, test_gc$statistic, test_gc$p.value, test_gc$alternative), c(test_depth$data.name, test_depth$method, test_depth$statistic, test_depth$p.value, test_depth$alternative), c(test_rep$data.name, test_rep$method, test_rep$statistic, test_rep$p.value, test_rep$alternative), c(test_genes$data.name, test_genes$method, test_genes$statistic, test_genes$p.value,test_genes$alternative), c(test_rs$data.name, test_rs$method, test_rs$statistic, test_rs$p.value, test_rs$alternative))))

		rownames(test_result_df) <- NULL
		test_result_df <- cbind(chr_name, test_result_df)
		colnames(test_result_df) <- c('chr_name', 'formula', 'test', 'statistic', 'p_value', 'alternative_hypothesis')
	}
	return(test_result_df)
}



all_chr_test_df <- rbind(
	cbind(test_significant_difference('chr 1', t.test), test_significant_difference('chr 1', wilcox.test)[,3:6], test_significant_difference('chr 1', var.test)[,3:6], test_significant_difference('chr 1', ks.test)[,3:6]),
	cbind(test_significant_difference('chr 2', t.test), test_significant_difference('chr 2', wilcox.test)[,3:6], test_significant_difference('chr 2', var.test)[,3:6], test_significant_difference('chr 1', ks.test)[,3:6]),
	cbind(test_significant_difference('PAR', t.test), test_significant_difference('PAR', wilcox.test)[,3:6], test_significant_difference('PAR', var.test)[,3:6], test_significant_difference('chr 1', ks.test)[,3:6])
)
write.table('all_statistical_comparisons.tsv', sep='\t', row.names=F, col.names=T, x=all_chr_test_df, quote=F)


####################################################################
# The histogram plot.
####################################################################
# This is the scientific way to determine plot limits, but it makes the x axis harder to read.
#centro_min_gc <- min(unlist(lapply(list(chr_1_short, chr_1_long, chr_2_short, chr_2_long, chr_PAR_short, chr_PAR_long), function(x) x$GC)))
#centro_max_gc <- max(unlist(lapply(list(chr_1_short, chr_1_long, chr_2_short, chr_2_long, chr_PAR_short, chr_PAR_long), function(x) x$GC)))
#centro_min_rep <- min(unlist(lapply(list(chr_1_short, chr_1_long, chr_2_short, chr_2_long, chr_PAR_short, chr_PAR_long), function(x) x$PctMasked)))
#centro_max_rep <- max(unlist(lapply(list(chr_1_short, chr_1_long, chr_2_short, chr_2_long, chr_PAR_short, chr_PAR_long), function(x) x$PctMasked)))
#centro_min_cov <- min(unlist(lapply(list(chr_1_short, chr_1_long, chr_2_short, chr_2_long, chr_PAR_short, chr_PAR_long), function(x) x$Depth)))
#centro_max_cov <- max(unlist(lapply(list(chr_1_short, chr_1_long, chr_2_short, chr_2_long, chr_PAR_short, chr_PAR_long), function(x) x$Depth)))
#centro_min_rs <- min(unlist(lapply(list(chr_1_short, chr_1_long, chr_2_short, chr_2_long, chr_PAR_short, chr_PAR_long), function(x) x$RESites)))
#centro_max_rs <- max(unlist(lapply(list(chr_1_short, chr_1_long, chr_2_short, chr_2_long, chr_PAR_short, chr_PAR_long), function(x) x$RESites)))

# This is the hack way to determine plot limits, which makes it look better. So who's the hack?
centro_min_gc <- 34
centro_max_gc <- 50
centro_min_rep <- 0
centro_max_rep <- 100
centro_min_cov <- 0
centro_max_cov <- 1200
centro_min_rs <- 0
centro_max_rs <- 350


pdf('centro_stats_gc_histograms.pdf')
par(mfrow=c(2, 3))
hist(chr_1_short$GC, xlab='GC', breaks=20, xlim=c(centro_min_gc, centro_max_gc), col='#e78ac3', main='Chr. 1 short arm')
abline(v=mean(chr_1_short$GC), col='red')
legend('topright', legend = list(paste('Mean:', round(mean(chr_1_short$GC), 2)), paste('SD:', round(sd(chr_1_short$GC), 2))), bty='n' )

hist(chr_2_short$GC, xlab='GC', breaks=20, xlim=c(centro_min_gc, centro_max_gc), col='#a6d854', main='Chr. 2 short arm')
abline(v=mean(chr_2_short$GC), col='red')
legend('topright', legend = list(paste('Mean:', round(mean(chr_2_short$GC), 2)), paste('SD:', round(sd(chr_PAR_short$GC), 2))), bty='n' )

hist(chr_PAR_short$GC, xlab='GC', breaks=20, xlim=c(centro_min_gc, centro_max_gc), col='#ffd92f', main='PAR short arm')
abline(v=mean(chr_PAR_short$GC), col='red')
legend('topright', legend = list(paste('Mean:', round(mean(chr_PAR_short$GC), 2)), paste('SD:', round(sd(chr_PAR_short$GC), 2))), bty='n' )

hist(chr_1_long$GC, xlab='GC', breaks=20, xlim=c(centro_min_gc, centro_max_gc), col='#f0b3d8', main='Chr.1 long arm')
abline(v=mean(chr_1_long$GC), col='red')
legend('topright', legend = list(paste('Mean:', round(mean(chr_1_long$GC), 2)), paste('SD:', round(sd(chr_1_long$GC), 2))), bty='n' )

hist(chr_2_long$GC, xlab='GC', breaks=20, xlim=c(centro_min_gc, centro_max_gc), col='#cce89e', main='Chr. 2 long arm')
abline(v=mean(chr_2_long$GC), col='red')
legend('topright', legend = list(paste('Mean:', round(mean(chr_2_long$GC), 2)), paste('SD:', round(sd(chr_2_long$GC), 2))), bty='n' )

hist(chr_PAR_long$GC, xlab='GC', breaks=20, xlim=c(centro_min_gc, centro_max_gc), col='#ffe883', main='PAR long arm')
abline(v=mean(chr_PAR_long$GC), col='red')
legend('topright', legend = list(paste('Mean:', round(mean(chr_PAR_long$GC), 2)), paste('SD:', round(sd(chr_PAR_long$GC), 2))), bty='n' )
dev.off()



####################################################################
# boxplots.
####################################################################
split_plot = FALSE
if(split_plot==TRUE){
	pdf('centro_stats_boxplots_split.pdf')
	par(mar=c(8,4.5,2,2))
	boxplot(x=list(chr_1_short$GC, chr_1_long$GC, chr_2_short$GC, chr_2_long$GC, chr_PAR_short$GC, chr_PAR_long$GC), pch=21, bg='white', col = c('#e78ac3', '#f0b3d8', '#a6d854', '#cce89e', '#ffd92f', '#ffe883'), names=c('Chr. 1 short arm', 'Chr. 1 long arm', 'Chr. 2 short arm', 'Chr. 2 long arm', 'PAR short arm', 'PAR long arm'), las=2, ylab='GC', main='GC content vs. chromosomal position', ylim=c(centro_min_gc, centro_max_gc))
	text(x=c(1.5, 3.5, 5.5), y=c(centro_max_gc-1, centro_max_gc-1, centro_max_gc-1), labels = c('*', '*', '*'))
	lines(x=c(1,2), y=c(centro_max_gc-1.5,centro_max_gc-1.5), lty=3)
	lines(x=c(3,4), y=c(centro_max_gc-1.5,centro_max_gc-1.5), lty=3)
	lines(x=c(5,6), y=c(centro_max_gc-1.5,centro_max_gc-1.5), lty=3)

	boxplot(x=list(chr_1_short$PctMasked, chr_1_long$PctMasked, chr_2_short$PctMasked, chr_2_long$PctMasked, chr_PAR_short$PctMasked, chr_PAR_long$PctMasked), pch=21, bg='white', col = c('#e78ac3', '#f0b3d8', '#a6d854', '#cce89e', '#ffd92f', '#ffe883'), names=c('Chr. 1 short arm', 'Chr. 1 long arm', 'Chr. 2 short arm', 'Chr. 2 long arm', 'PAR short arm', 'PAR long arm'), las=2, ylab='Percent masked', main='Repetitive content vs. chromosomal position', ylim=c(centro_min_rep,centro_max_rep))
	text(x=c(1.5, 3.5, 5.5), y=c(centro_max_rep-2,centro_max_rep-2,centro_max_rep-2), labels = c('*', '*', '*'))
	lines(x=c(1,2), y=c(centro_max_rep-4,centro_max_rep-4), lty=3)
	lines(x=c(3,4), y=c(centro_max_rep-4,centro_max_rep-4), lty=3)
	lines(x=c(5,6), y=c(centro_max_rep-4,centro_max_rep-4), lty=3)

	boxplot(x=list(chr_1_short$gene_counts, chr_1_long$gene_counts, chr_2_short$gene_counts, chr_2_long$gene_counts, chr_PAR_short$gene_counts, chr_PAR_long$gene_counts), pch=21, bg='white', col = c('#e78ac3', '#f0b3d8', '#a6d854', '#cce89e', '#ffd92f', '#ffe883'), names=c('Chr. 1 short arm', 'Chr. 1 long arm', 'Chr. 2 short arm', 'Chr. 2 long arm', 'PAR short arm', 'PAR long arm'), las=2, ylab='Number of genes per 50kb', main='Gene density vs. chromosomal position')

	boxplot(x=list(chr_1_short$RESites, chr_1_long$RESites, chr_2_short$RESites, chr_2_long$RESites, chr_PAR_short$RESites, chr_PAR_long$RESites), pch=21, bg='white', col = c('#e78ac3', '#f0b3d8', '#a6d854', '#cce89e', '#ffd92f', '#ffe883'), names=c('Chr. 1 short arm', 'Chr. 1 long arm', 'Chr. 2 short arm', 'Chr. 2 long arm', 'PAR short arm', 'PAR long arm'), las=2, ylab='DpnII sites per 50kb', main='DpnII sites vs. chromosomal position', ylim=c(centro_min_rs,centro_max_rs))
	text(x=c(1.5, 3.5, 5.5), y=c(centro_max_rs-2,centro_max_rs-2,centro_max_rs-2), labels = c('*', '*', '*'))
	lines(x=c(1,2), y=c(centro_max_rs-4,centro_max_rs-4), lty=3)
	lines(x=c(3,4), y=c(centro_max_rs-4,centro_max_rs-4), lty=3)
	lines(x=c(5,6), y=c(centro_max_rs-4,centro_max_rs-4), lty=3)
	dev.off()
}
if(split_plot==FALSE){
	pdf('centro_stats_boxplots_joined_2.pdf')
	par(mar=c(1,6,2,2), mfrow=c(6,1))
	boxplot(x=list(chr_1_short$GC, chr_1_long$GC, chr_2_short$GC, chr_2_long$GC, chr_PAR_short$GC, chr_PAR_long$GC), pch=21, bg='white', col = c('#e78ac3', '#f0b3d8', '#a6d854', '#cce89e', '#ffd92f', '#ffe883'), las=2, ylab='GC\n(per 50kb)', ylim=c(centro_min_gc, centro_max_gc), xaxt='n', cex=1.5)
	text(x=c(1.5, 3.5, 5.5), y=c(centro_max_gc-1, centro_max_gc-1, centro_max_gc-1), labels = c('*', '*', '*'))
	lines(x=c(1,2), y=c(centro_max_gc-1.5,centro_max_gc-1.5), lty=3)
	lines(x=c(3,4), y=c(centro_max_gc-1.5,centro_max_gc-1.5), lty=3)
	lines(x=c(5,6), y=c(centro_max_gc-1.5,centro_max_gc-1.5), lty=3)

	boxplot(x=list(chr_1_short$Depth, chr_1_long$Depth, chr_2_short$Depth, chr_2_long$Depth, chr_PAR_short$Depth, chr_PAR_long$Depth), pch=21, bg='white', col = c('#e78ac3', '#f0b3d8', '#a6d854', '#cce89e', '#ffd92f', '#ffe883'), las=2, ylab='Coverage\n(per 50kb)', cex=1.5, xaxt='n')
	text(x=c(1.5), y=c(centro_max_cov), labels = c('*'))
	lines(x=c(1,2), y=c(centro_max_cov-4,centro_max_cov-4), lty=3)

	boxplot(x=list(chr_1_short$PctMasked, chr_1_long$PctMasked, chr_2_short$PctMasked, chr_2_long$PctMasked, chr_PAR_short$PctMasked, chr_PAR_long$PctMasked), pch=21, bg='white', col = c('#e78ac3', '#f0b3d8', '#a6d854', '#cce89e', '#ffd92f', '#ffe883'), las=2, ylab='Percent repetitive\n(per 50kb)', ylim=c(centro_min_rep,centro_max_rep), xaxt='n', cex=1.5)
	text(x=c(1.5, 3.5, 5.5), y=c(centro_max_rep-2,centro_max_rep-2,centro_max_rep-2), labels = c('*', '*', '*'))
	lines(x=c(1,2), y=c(centro_max_rep-4,centro_max_rep-4), lty=3)
	lines(x=c(3,4), y=c(centro_max_rep-4,centro_max_rep-4), lty=3)
	lines(x=c(5,6), y=c(centro_max_rep-4,centro_max_rep-4), lty=3)


	#par(mar=c(2,4.5,2,2))
	boxplot(x=list(chr_1_short$gene_counts, chr_1_long$gene_counts, chr_2_short$gene_counts, chr_2_long$gene_counts, chr_PAR_short$gene_counts, chr_PAR_long$gene_counts), pch=21, bg='white', col = c('#e78ac3', '#f0b3d8', '#a6d854', '#cce89e', '#ffd92f', '#ffe883'), las=2, ylab='Number of genes\n(per 50kb)', cex=1.5, xaxt='n')

	# DpnII sites
	boxplot(x=list(chr_1_short$RESites, chr_1_long$RESites, chr_2_short$RESites, chr_2_long$RESites, chr_PAR_short$RESites, chr_PAR_long$RESites), pch=21, bg='white', col = c('#e78ac3', '#f0b3d8', '#a6d854', '#cce89e', '#ffd92f', '#ffe883'), las=2, ylim=c(centro_min_rs, centro_max_rs), ylab='DpnII sites\n(per 50kb)', cex=1.5, names=c('Chr. 1 short arm', 'Chr. 1 long arm', 'Chr. 2 short arm', 'Chr. 2 long arm', 'PAR short arm', 'PAR long arm'))
	text(x=c(1.5, 3.5, 5.5), y=c(centro_max_rs-2,centro_max_rs-2,centro_max_rs-2), labels = c('*', '*', '*'))
	lines(x=c(1,2), y=c(centro_max_rs-4,centro_max_rs-4), lty=3)
	lines(x=c(3,4), y=c(centro_max_rs-4,centro_max_rs-4), lty=3)
	lines(x=c(5,6), y=c(centro_max_rs-4,centro_max_rs-4), lty=3)

	dev.off()
}



####################################################################
# Giant GC correlation plots. Probably not for publication.
####################################################################
# Chr. 1
pdf('centro_correlations.pdf')
par(mfrow=c(3,2), mar=c(4,5,4,2))
plot(chr_1_short$PctMasked, chr_1_short$GC,  xlab='Percent repetitive', ylab='GC', pch=21, bg=adjustcolor('#e78ac3', 0.3), col=NA, main='Chr. 1 short arm')
a <- lm(chr_1_short$GC ~ chr_1_short$PctMasked)
abline(a, lty=3)
ar <- summary(a)$adj.r.squared
ap <- summary(a)$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), list(MYVALUE = format(ar, dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), list(MYOTHERVALUE = format(ap, digits = 2)))[2]
legend('topright', legend=rp, bty='n')

plot(chr_1_long$PctMasked, chr_1_long$GC,  xlab='Percent repetitive', ylab='GC', pch=21, bg=adjustcolor('#f0b3d8', 0.3), col=NA, main='Chr. 1 long arm')
a <- lm(chr_1_long$GC ~ chr_1_long$PctMasked)
abline(a, lty=3)
ar <- summary(a)$adj.r.squared
ap <- summary(a)$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), list(MYVALUE = format(ar, dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), list(MYOTHERVALUE = format(ap, digits = 2)))[2]
legend('topright', legend=rp, bty='n')

plot(chr_1_short$Depth, chr_1_short$GC,  xlab='Depth', ylab='GC', pch=21, bg=adjustcolor('#e78ac3', 0.3), col=NA)
a <- lm(chr_1_short$GC ~ chr_1_short$Depth)
abline(a, lty=3)
ar <- summary(a)$adj.r.squared
ap <- summary(a)$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), list(MYVALUE = format(ar, dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), list(MYOTHERVALUE = format(ap, digits = 2)))[2]
legend('topright', legend=rp, bty='n')

plot(chr_1_long$Depth, chr_1_long$GC,  xlab='Depth', ylab='GC', pch=21, bg=adjustcolor('#f0b3d8', 0.3), col=NA)
a <- lm(chr_1_long$GC ~ chr_1_long$Depth)
abline(a, lty=3)
ar <- summary(a)$adj.r.squared
ap <- summary(a)$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), list(MYVALUE = format(ar, dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), list(MYOTHERVALUE = format(ap, digits = 2)))[2]
legend('topright', legend=rp, bty='n')

plot(chr_1_short$gene_counts, chr_1_short$GC,  xlab='Genes per 5kb', ylab='GC', pch=21, bg=adjustcolor('#e78ac3', 0.3), col=NA)
a <- lm(chr_1_short$GC ~ chr_1_short$gene_counts)
abline(a, lty=3)
ar <- summary(a)$adj.r.squared
ap <- summary(a)$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), list(MYVALUE = format(ar, dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), list(MYOTHERVALUE = format(ap, digits = 2)))[2]
legend('topright', legend=rp, bty='n')

plot(chr_1_long$gene_counts, chr_1_long$GC,  xlab='Genes per 5kb', ylab='GC', pch=21, bg=adjustcolor('#f0b3d8', 0.3), col=NA)
a <- lm(chr_1_long$GC ~ chr_1_long$gene_counts)
abline(a, lty=3)
ar <- summary(a)$adj.r.squared
ap <- summary(a)$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), list(MYVALUE = format(ar, dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), list(MYOTHERVALUE = format(ap, digits = 2)))[2]
legend('topright', legend=rp, bty='n')

# Chr. 2
plot(chr_2_short$PctMasked, chr_2_short$GC,  xlab='Percent repetitive', ylab='GC', pch=21, bg=adjustcolor('#a6d854', 0.3), col=NA, main='Chr. 2 short arm')
a <- lm(chr_2_short$GC ~ chr_2_short$PctMasked)
abline(a, lty=3)
ar <- summary(a)$adj.r.squared
ap <- summary(a)$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), list(MYVALUE = format(ar, dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), list(MYOTHERVALUE = format(ap, digits = 2)))[2]
legend('topright', legend=rp, bty='n')

plot(chr_2_long$PctMasked, chr_2_long$GC,  xlab='Percent repetitive', ylab='GC', pch=21, bg=adjustcolor('#cce89e', 0.3), col=NA, main='Chr. 2 long arm')
a <- lm(chr_2_long$GC ~ chr_2_long$PctMasked)
abline(a, lty=3)
ar <- summary(a)$adj.r.squared
ap <- summary(a)$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), list(MYVALUE = format(ar, dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), list(MYOTHERVALUE = format(ap, digits = 2)))[2]
legend('topright', legend=rp, bty='n')

plot(chr_2_short$Depth, chr_2_short$GC,  xlab='Depth', ylab='GC', pch=21, bg=adjustcolor('#a6d854', 0.3), col=NA)
a <- lm(chr_2_short$GC ~ chr_2_short$Depth)
abline(a, lty=3)
ar <- summary(a)$adj.r.squared
ap <- summary(a)$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), list(MYVALUE = format(ar, dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), list(MYOTHERVALUE = format(ap, digits = 2)))[2]
legend('topright', legend=rp, bty='n')

plot(chr_2_long$Depth, chr_2_long$GC,  xlab='Depth', ylab='GC', pch=21, bg=adjustcolor('#cce89e', 0.3), col=NA)
a <- lm(chr_2_long$GC ~ chr_2_long$Depth)
abline(a, lty=3)
ar <- summary(a)$adj.r.squared
ap <- summary(a)$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), list(MYVALUE = format(ar, dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), list(MYOTHERVALUE = format(ap, digits = 2)))[2]
legend('topright', legend=rp, bty='n')

plot(chr_2_short$gene_counts, chr_2_short$GC,  xlab='Genes per 5kb', ylab='GC', pch=21, bg=adjustcolor('#a6d854', 0.3), col=NA)
a <- lm(chr_2_short$GC ~ chr_2_short$gene_counts)
abline(a, lty=3)
ar <- summary(a)$adj.r.squared
ap <- summary(a)$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), list(MYVALUE = format(ar, dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), list(MYOTHERVALUE = format(ap, digits = 2)))[2]
legend('topright', legend=rp, bty='n')

plot(chr_2_long$gene_counts, chr_2_long$GC,  xlab='Genes per 5kb', ylab='GC', pch=21, bg=adjustcolor('#cce89e', 0.3), col=NA)
a <- lm(chr_2_long$GC ~ chr_2_long$gene_counts)
abline(a, lty=3)
ar <- summary(a)$adj.r.squared
ap <- summary(a)$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), list(MYVALUE = format(ar, dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), list(MYOTHERVALUE = format(ap, digits = 2)))[2]
legend('topright', legend=rp, bty='n')

# PAR
plot(chr_PAR_short$PctMasked, chr_PAR_short$GC,  xlab='Percent repetitive', ylab='GC', pch=21, bg=adjustcolor('#ffd92f', 0.3), col=NA, main='PAR short arm')
a <- lm(chr_PAR_short$GC ~ chr_PAR_short$PctMasked)
abline(a, lty=3)
ar <- summary(a)$adj.r.squared
ap <- summary(a)$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), list(MYVALUE = format(ar, dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), list(MYOTHERVALUE = format(ap, digits = 2)))[2]
legend('topright', legend=rp, bty='n')

plot(chr_PAR_long$PctMasked, chr_PAR_long$GC,  xlab='Percent repetitive', ylab='GC', pch=21, bg=adjustcolor('#ffe883', 0.3), col=NA, main='PAR long arm')
a <- lm(chr_PAR_long$GC ~ chr_PAR_long$PctMasked)
abline(a, lty=3)
ar <- summary(a)$adj.r.squared
ap <- summary(a)$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), list(MYVALUE = format(ar, dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), list(MYOTHERVALUE = format(ap, digits = 2)))[2]
legend('topright', legend=rp, bty='n')

plot(chr_PAR_short$Depth, chr_PAR_short$GC,  xlab='Depth', ylab='GC', pch=21, bg=adjustcolor('#ffd92f', 0.3), col=NA)
a <- lm(chr_PAR_short$GC ~ chr_PAR_short$Depth)
abline(a, lty=3)
ar <- summary(a)$adj.r.squared
ap <- summary(a)$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), list(MYVALUE = format(ar, dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), list(MYOTHERVALUE = format(ap, digits = 2)))[2]
legend('topright', legend=rp, bty='n')

plot(chr_PAR_long$Depth, chr_PAR_long$GC,  xlab='Depth', ylab='GC', pch=21, bg=adjustcolor('#ffe883', 0.3), col=NA)
a <- lm(chr_PAR_long$GC ~ chr_PAR_long$Depth)
abline(a, lty=3)
ar <- summary(a)$adj.r.squared
ap <- summary(a)$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), list(MYVALUE = format(ar, dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), list(MYOTHERVALUE = format(ap, digits = 2)))[2]
legend('topright', legend=rp, bty='n')

plot(chr_PAR_short$gene_counts, chr_PAR_short$GC,  xlab='Genes per 5kb', ylab='GC', pch=21, bg=adjustcolor('#ffd92f', 0.3), col=NA)
a <- lm(chr_PAR_short$GC ~ chr_PAR_short$gene_counts)
abline(a, lty=3)
ar <- summary(a)$adj.r.squared
ap <- summary(a)$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), list(MYVALUE = format(ar, dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), list(MYOTHERVALUE = format(ap, digits = 2)))[2]
legend('topright', legend=rp, bty='n')

plot(chr_PAR_long$gene_counts, chr_PAR_long$GC,  xlab='Genes per 5kb', ylab='GC', pch=21, bg=adjustcolor('#ffe883', 0.3), col=NA)
a <- lm(chr_PAR_long$GC ~ chr_PAR_long$gene_counts)
abline(a, lty=3)
ar <- summary(a)$adj.r.squared
ap <- summary(a)$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), list(MYVALUE = format(ar, dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), list(MYOTHERVALUE = format(ap, digits = 2)))[2]
legend('topright', legend=rp, bty='n')

dev.off()
