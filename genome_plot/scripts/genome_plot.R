#!/bin/Rscript

library(data.table)
library(parallel)
library(optparse)

####################################################################
#
# Options for inputs and outputs. Grabs a few things so that the
# script can point to the right files and generate the requested
# plots.
#
####################################################################

option_list <- list(
  make_option(c('-F', '--windowFile'), type='character', default=NA, help='String. Path to file containing windowed sequence statistics (i.e., the "windowed_stats.txt" output of genome_statistics.py). No default.'),
	make_option(c('-G', '--geneFile'), type='character', default=NA, help='String. GFF file containing gene coordinates. No default.'),
	make_option(c('-A', '--gapFile'), type='character', default=NA, help='String. GFF file containing gap coordinates. No default.'),
	make_option(c('-O', '--outputFile'), type='character', default=NA, help='String. Path to output PDF file. No default.'),
	make_option(c('-W', '--windowSize'), type='integer', default=NA, help='Integer. Size of the window used to generate the windowed_stats.txt file. No default.'),
	make_option(c('-s', '--smoothSplines'), action='store_true', default=FALSE, help='Boolean. Use smoothing splines in visualization. Default: FALSE'),
	make_option(c('-k', '--scaledPlot'), action='store_true', default=FALSE, help='Boolean. Make each track scaled relative to each chromosome. Default: FALSE'),
	make_option(c('-m', '--minimumZero'), action='store_true', default=FALSE, help='Boolean. Do not allow splines to go below 0. Default: FALSE'),
  make_option(c('-t', '--threads'), type='integer', default=1, help='Integer. Number of threads. Default: 1')
)

parser <- OptionParser(usage="%prog [-F windowFile] [-G geneFile] [-A gapfile] [-o outputFile] [-W windowSize] [additional optional arguments]", option_list=option_list)

args <- parse_args(parser)

# Gather input options
window_file <- args$windowFile
gene_file <- args$geneFile
gap_file <- args$gapFile
output_file <- args$outputFile
window_size <- args$windowSize

# Make sure that all of the required input options actually exist, and
# throw an error and stop the script if they are missing.
check_required_options <- sapply(names(args), simplify=FALSE, USE.NAMES=TRUE, function(arg_name) {
	arg_value <- args[[arg_name]]
	if(is.na(arg_value)){
		stop(paste('Error: required input --', arg_name, ' is not specified.', sep=''))
	}
})

# If the argument for output_file does not have pdf file extension, add it.
output_file_ext <- tail(strsplit(output_file, split='\\.')[[1]], n=1)
if(output_file_ext != 'pdf'){
	output_file <- paste(output_file, '.pdf', sep='')
}

# Check all of the optional inputs and set their defaults if the defaults do not yet exist.
# Smooth splines?
if(isTRUE(args$smoothSplines)){
	plot_smooth = TRUE
} else {
	plot_smooth = FALSE
}

# Scaled tracks per-chromosome?
if(isTRUE(args$scaledPlot)){
	plot_scaled = TRUE
} else {
	plot_scaled = FALSE
}

# Allow splines below 0?
if(isTRUE(args$minimumZero)){
	plot_replace_zeroes = TRUE
} else {
	plot_replace_zeroes = FALSE
}

# Number of threads?
num_threads <- args$threads
if(num_threads > 1){
	make_cluster <- TRUE
	setDTthreads(threads = num_threads)
} else {
	make_cluster <- FALSE
	setDTthreads(1)
}

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

# Master plotting function. Complicated to allow for different parameter combinations.
unified_chr_plot <- function(scaffold_name, smoothed_splines=FALSE, scaled_per_chromosome=FALSE, overwrite_zeroes=FALSE){
	# Get data for a given chromosome
  gene_data <- split_genes[scaffold_name][[1]]
  scaff_data <- split_stats[scaffold_name][[1]]
  gap_data <- split_gaps[scaffold_name][[1]]
  chr_size <- scaffold_sizes$ends[match(scaffold_name, scaffold_sizes$chr_name)]

	######################################################
  # Making the plot
  ######################################################
	# Ranges from 0 to the max chromosome size on X, and from the bottom to the top of the individual chromsome's fitted spline
	if(scaled_per_chromosome==TRUE && smoothed_splines==TRUE){
	  # GC content plot
		plot(x=c(0, max(scaffold_sizes$ends)), y=c(min(spline_coords(scaff_data$Window_Mid,scaff_data$GC,makeNonZero=overwrite_zeroes)), max(spline_coords(scaff_data$Window_Mid,scaff_data$GC,makeNonZero=overwrite_zeroes))), xlab=NA, ylab='GC', type='n', las=2, xaxt='none')
	  polygon(x=c(0, scaff_data$Window_Mid, chr_size), y=c(min(spline_coords(scaff_data$Window_Mid,scaff_data$GC,makeNonZero=overwrite_zeroes)), spline_coords(scaff_data$Window_Mid,scaff_data$GC,makeNonZero=overwrite_zeroes), min(spline_coords(scaff_data$Window_Mid,scaff_data$GC,makeNonZero=overwrite_zeroes))), col='#fdbe85', border='#fd8d3c')

		# Coverage plot
		# As above, x ranges from - to max chrosomome. Y ranges from bottom to top of spline
		plot(x=c(0, max(scaffold_sizes$ends)), y=c(min(spline_coords(scaff_data$Window_Mid,scaff_data$Depth,makeNonZero=overwrite_zeroes)), max(spline_coords(scaff_data$Window_Mid,scaff_data$Depth,makeNonZero=overwrite_zeroes))), xlab=NA, ylab='Coverage', type='n', las=2, xaxt='none')
	  polygon(x=c(0, scaff_data$Window_Mid, chr_size), y=c(min(spline_coords(scaff_data$Window_Mid,scaff_data$Depth,makeNonZero=overwrite_zeroes)), spline_coords(scaff_data$Window_Mid,scaff_data$Depth,makeNonZero=overwrite_zeroes), min(spline_coords(scaff_data$Window_Mid,scaff_data$Depth,makeNonZero=overwrite_zeroes))), col='#bdd7e7', border='#6baed6')

		# Repetitive region plot
	  plot(x=c(0, max(scaffold_sizes$ends)), y=c(min(spline_coords(scaff_data$Window_Mid,scaff_data$PctMasked,makeNonZero=overwrite_zeroes)), max(spline_coords(scaff_data$Window_Mid,scaff_data$PctMasked,makeNonZero=overwrite_zeroes))), xlab=NA, ylab='Percent repetitive', type='n', las=2, xaxt='none')
	  polygon(x=c(0,scaff_data$Window_Mid,chr_size), y=c(min(spline_coords(scaff_data$Window_Mid,scaff_data$PctMasked,makeNonZero=overwrite_zeroes)), spline_coords(scaff_data$Window_Mid,scaff_data$PctMasked,makeNonZero=overwrite_zeroes), min(spline_coords(scaff_data$Window_Mid,scaff_data$PctMasked,makeNonZero=overwrite_zeroes))), col='#cbc9e2', border='#9e9ac8')

		# Gene density plot
	  plot(x=c(0, max(scaffold_sizes$ends)), y=c(min(spline_coords(scaff_data$Window_Mid,scaff_data$gene_counts,makeNonZero=overwrite_zeroes)), max(spline_coords(scaff_data$Window_Mid,scaff_data$gene_counts,makeNonZero=overwrite_zeroes))), xlab=scaffold_name, ylab='Gene density', type='n', las=2, xaxt='none')
	  axis(1, at=c(0, chr_size*0.33, chr_size*0.66, chr_size), labels=paste(c(0, round(chr_size/1000000*0.33, 2), round(chr_size/1000000*0.66, 2), round(chr_size/1000000, 2)), "Mb"), las = 2)
	  # Annotate gap regions
	  polygon(x=c(0, scaff_data$Window_Mid, chr_size), y=c(min(spline_coords(scaff_data$Window_Mid,scaff_data$gene_counts,makeNonZero=overwrite_zeroes)), spline_coords(scaff_data$Window_Mid,scaff_data$gene_counts,makeNonZero=overwrite_zeroes), min(spline_coords(scaff_data$Window_Mid,scaff_data$gene_counts,makeNonZero=overwrite_zeroes))), col='#ffff99', border='#FFE88C')
		# Put gap rectangles on top of the polygon
		rect(xleft=gap_data$start, xright=gap_data$end, ybottom=rep(min(spline_coords(scaff_data$Window_Mid,scaff_data$gene_counts,makeNonZero=overwrite_zeroes)), nrow(gap_data)), ytop=rep(max(spline_coords(scaff_data$Window_Mid,scaff_data$gene_counts,makeNonZero=overwrite_zeroes)), nrow(gap_data)), col='#fb8072', border='#fb8072')

	} else if(scaled_per_chromosome==TRUE && smoothed_splines==FALSE) {
		# GC content plot
		plot(x=c(0, max(scaffold_sizes$ends)), y=c(min(scaff_data$GC), max(scaff_data$GC)), xlab=NA, ylab='GC', type='n', las=2, xaxt='none')
		polygon(x=c(0, scaff_data$Window_Mid, chr_size), y=c(min(scaff_data$GC), scaff_data$GC, min(scaff_data$GC)), col='#fdbe85', border='#fd8d3c')

		# Coverage plot
		plot(x=c(0, max(scaffold_sizes$ends)), y=c(min(scaff_data$Depth), max(scaff_data$Depth)), xlab=NA, ylab='Coverage', type='n', las=2, xaxt='none')
		polygon(x=c(0, scaff_data$Window_Mid, chr_size), y=c(min(scaff_data$Depth), scaff_data$Depth, min(scaff_data$Depth)), col='#bdd7e7', border='#6baed6')

		# Repetitive region plot
	  plot(x=c(0, max(scaffold_sizes$ends)), y=c(min(scaff_data$PctMasked), max(scaff_data$PctMasked)), xlab=NA, ylab='Percent repetitive', type='n', las=2, xaxt='none')
	  polygon(x=c(0,scaff_data$Window_Mid,chr_size), y=c(min(scaff_data$PctMasked), scaff_data$PctMasked, min(scaff_data$PctMasked)), col='#cbc9e2', border='#9e9ac8')

		# Gene density plot
	  plot(x=c(0, max(scaffold_sizes$ends)), y=c(min(scaff_data$gene_counts), max(scaff_data$gene_counts)), xlab=scaffold_name, ylab='Gene density', type='n', las=2, xaxt='none')
	  axis(1, at=c(0, chr_size*0.33, chr_size*0.66, chr_size), labels=paste(c(0, round(chr_size/1000000*0.33, 2), round(chr_size/1000000*0.66, 2), round(chr_size/1000000, 2)), "Mb"), las = 2)
	  # Annotate gap regions
	  polygon(x=c(0, scaff_data$Window_Mid, chr_size), y=c(min(scaff_data$gene_counts), scaff_data$gene_counts, min(scaff_data$gene_counts)), col='#ffff99', border='#FFE88C')
		# Put gap rectangles on top of the polygon
		rect(xleft=gap_data$start, xright=gap_data$end, ybottom=rep(min(scaff_data$gene_counts), nrow(gap_data)), ytop=rep(max(scaff_data$gene_counts), nrow(gap_data)), col='#fb8072', border='#fb8072')

	} else if (scaled_per_chromosome==FALSE && smoothed_splines==TRUE){
		# GC content plot
		plot(x=c(0, max(scaffold_sizes$ends)), y=c(overall_gc_min, overall_gc_max), xlab=NA, ylab='GC', type='n', las=2, xaxt='none')
	  polygon(x=c(0, scaff_data$Window_Mid, chr_size), y=c(overall_gc_min, spline_coords(scaff_data$Window_Mid,scaff_data$GC,makeNonZero=overwrite_zeroes), overall_gc_min), col='#fdbe85', border='#fd8d3c')

		# Coverage plot
		plot(x=c(0, max(scaffold_sizes$ends)), y=c(overall_cov_min, overall_cov_max), xlab=NA, ylab='Coverage', type='n', las=2, xaxt='none')
	  polygon(x=c(0, scaff_data$Window_Mid, chr_size), y=c(overall_cov_min, spline_coords(scaff_data$Window_Mid,scaff_data$Depth,makeNonZero=overwrite_zeroes), overall_cov_min), col='#bdd7e7', border='#6baed6')

		# Repetitive region plot
	  plot(x=c(0, max(scaffold_sizes$ends)), y=c(overall_rep_min, overall_rep_max), xlab=NA, ylab='Percent repetitive', type='n', las=2, xaxt='none')
	  polygon(x=c(0,scaff_data$Window_Mid,chr_size), y=c(overall_rep_min, spline_coords(scaff_data$Window_Mid,scaff_data$PctMasked,makeNonZero=overwrite_zeroes), overall_rep_min), col='#cbc9e2', border='#9e9ac8')

		# Gene density plot
	  plot(x=c(0, max(scaffold_sizes$ends)), y=c(overall_gene_min, overall_gene_max), xlab=scaffold_name, ylab='Gene density', type='n', las=2, xaxt='none')
	  axis(1, at=c(0, chr_size*0.33, chr_size*0.66, chr_size), labels=paste(c(0, round(chr_size/1000000*0.33, 2), round(chr_size/1000000*0.66, 2), round(chr_size/1000000, 2)), "Mb"), las = 2)
	  # Annotate gap regions
	  polygon(x=c(0, scaff_data$Window_Mid, chr_size), y=c(overall_gene_min, spline_coords(scaff_data$Window_Mid,scaff_data$gene_counts,makeNonZero=overwrite_zeroes), overall_gene_min), col='#ffff99', border='#FFE88C')
		# Put gap rectangles on top of the polygon
		rect(xleft=gap_data$start, xright=gap_data$end, ybottom=rep(overall_gene_min, nrow(gap_data)), ytop=rep(overall_gene_max, nrow(gap_data)), col='#fb8072', border='#fb8072')
	} else if (scaled_per_chromosome==FALSE && smoothed_splines==FALSE){
		# GC content plot
		plot(x=c(0, max(scaffold_sizes$ends)), y=c(min(do.call(rbind, split_stats)$GC), max(do.call(rbind, split_stats)$GC)), xlab=NA, ylab='GC', type='n', las=2, xaxt='none')
		polygon(x=c(0, scaff_data$Window_Mid, chr_size), y=c(min(do.call(rbind, split_stats)$GC), scaff_data$GC, min(do.call(rbind, split_stats)$GC)), col='#fdbe85', border='#fd8d3c')

		# Coverage plot
		plot(x=c(0, max(scaffold_sizes$ends)), y=c(min(do.call(rbind, split_stats)$Depth), max(do.call(rbind, split_stats)$Depth)), xlab=NA, ylab='Coverage', type='n', las=2, xaxt='none')
		polygon(x=c(0, scaff_data$Window_Mid, chr_size), y=c(min(do.call(rbind, split_stats)$Depth), scaff_data$Depth, min(do.call(rbind, split_stats)$Depth)), col='#bdd7e7', border='#6baed6')

		# Repetitive region plot
		plot(x=c(0, max(scaffold_sizes$ends)), y=c(min(do.call(rbind, split_stats)$PctMasked), max(do.call(rbind, split_stats)$PctMasked)), xlab=NA, ylab='Percent repetitive', type='n', las=2, xaxt='none')
	  polygon(x=c(0,scaff_data$Window_Mid,chr_size), y=c(min(do.call(rbind, split_stats)$PctMasked), scaff_data$PctMasked, min(do.call(rbind, split_stats)$PctMasked)), col='#cbc9e2', border='#9e9ac8')

		# Gene density plot
		plot(x=c(0, max(scaffold_sizes$ends)), y=c(min(do.call(rbind, split_stats)$gene_counts), max(do.call(rbind, split_stats)$gene_counts)), xlab=scaffold_name, ylab='Gene density', type='n', las=2, xaxt='none')
	  axis(1, at=c(0, chr_size*0.33, chr_size*0.66, chr_size), labels=paste(c(0, round(chr_size/1000000*0.33, 2), round(chr_size/1000000*0.66, 2), round(chr_size/1000000, 2)), "Mb"), las = 2)
	  # Annotate gap regions
	  polygon(x=c(0, scaff_data$Window_Mid, chr_size), y=c(min(do.call(rbind, split_stats)$gene_counts), scaff_data$gene_counts, min(do.call(rbind, split_stats)$gene_counts)), col='#ffff99', border='#FFE88C')
		rect(xleft=gap_data$start, xright=gap_data$end, ybottom=rep(min(do.call(rbind, split_stats)$gene_counts), nrow(gap_data)), ytop=rep(max(do.call(rbind, split_stats)$gene_counts), nrow(gap_data)), col='#fb8072', border='#fb8072')
	}
}

####################################################################
#
# The plotting steps. The type of plot is specified by the input
# options. The plotting parameters are complicated, but the plot
# is pretty complicated, so...
#
####################################################################

pdf(output_file)
# These plotting options are always the same for every permutation of the plots.
# This is how we set track sizes, margin sizes, etc.
par(mar=c(1.2,4,1,1), cex.axis=0.4, cex.lab=0.4, mgp=c(3,0.5,0))
layout(matrix(c(1,1,5,5, 2,2,6,6, 3,3,7,7, 4,4,8,8, 9,9,13,13, 10,10,14,14, 11,11,15,15, 12,12,16,16, 17,17,21,21, 18,18,22,22, 19,19,23,23, 20,20,24,24), nrow=12, byrow=T))
# Now, determine the plotting options, and make the plot for every chromosome.
lapply(list('chr 1', 'chr 2', 'PAR', 'XSR', 'YSR'), function(x) {
	unified_chr_plot(scaffold_name=x, smoothed_splines=plot_smooth, scaled_per_chromosome=plot_scaled, overwrite_zeroes=plot_replace_zeroes)
} )
dev.off()
