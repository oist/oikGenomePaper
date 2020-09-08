intergenic_spaces <- function(genes_gff){
	intergenic_spaces <- unlist(lapply(1:nrow(genes_gff), function(X){
		gene_end <- as.numeric(genes_gff[X,5])
		next_gene <- as.numeric(genes_gff[X+1, 4])

		if(!is.na(next_gene)){
			difference <- next_gene - gene_end
			if(difference >= 0){
				return(difference)
			}
		}
	}))
	return(intergenic_spaces)
}

intergenic_spaces(split_genes[[1]])

