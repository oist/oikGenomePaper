library(readxl)
library(tidyverse)
library(VennDiagram)
library(eulerr)

venn_diagram <- read_excel("~/Dropbox (OIST)/Aki's Documents/R/BUSCO/BUSCO_Missing_IDs_per_genome.xlsx")

Oki <- venn_diagram$`I69-3`
Oki <- venn_diagram$`I69-3` %>% discard(is.na)
Nor <- venn_diagram$Norway  
Osa <- venn_diagram$OSKA2016 %>% discard(is.na)  
l <- list(Oki=Oki, Osa=Osa, Nor=Nor)
sapply(l, length)
venn <- get.venn.partitions(l, keep.elements = FALSE)
venn
plot(euler(l), shape = "ellipse", quantities = TRUE)
