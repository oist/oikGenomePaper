library(VennDiagram)
library(readxl)
library(tidyverse)
library(RColorBrewer)
library(data.table)
library(magrittr)
library(reshape2)

venn_diagram <- read_excel("~/Dropbox (OIST)/Aki's Documents/R/BUSCO/BUSCO_Missing_IDs_per_genome.xlsx")

Oki <-venn_diagram[ ,1] 
Nor <- venn_diagram[ ,2]
Osa <- venn_diagram[ ,3]

Pastel2 <- brewer.pal(3, "Pastel2")
Set2 <- brewer.pal(3, "Set2")
Accent <- brewer.pal(3, "Accent")

venn.diagram(
  x = list(Oki, Osa, Nor),
  category.names = c("Okinawa" , "Osaka" , "Norway"),
  filename = 'BUSCO_venn_diagram.tiff',
  imagetype = "tiff",
  na = "none",
  cex = 1,
  fontfamily = "sans",
  cat.cex = 1.5,
  cat.fontfamily = "sans",
  col = c("#440154ff", '#21908dff', '#fde725ff'),
  cat.col = c("#440154ff", '#21908dff', '#fde725ff'),
  output=TRUE
)


# count the number of occurrence 
no_NA <- venn_diagram %>%
  gather(variable, value) %>%
  drop_na() 

count <- no_NA %>%
  add_count(variable, name = "total # of genes") %>%
  add_count(value, sort = T, name = "matched") %>%
  group_by(variable) %>%
  add_count(matched, name = "# of matched genes") 

match_all <- dcast(no_NA, value ~ variable, fun.aggregate = length)

all_3 <- match_all %>% filter(`I69-3` == "1" & Norway == "1" & OSKA2016 == "1") %>% add_count
Oki_Nor <- match_all %>% filter(`I69-3` =="1" & Norway =="1" & OSKA2016 == "0") %>% add_count()
Oki_Osa <- match_all %>% filter(`I69-3` =="1" & Norway =="0" & OSKA2016 == "1") %>% add_count()
Nor_Osa <- match_all %>% filter(`I69-3` =="0" & Norway =="1" & OSKA2016 == "1") %>% add_count()
Oki_only <- match_all %>% filter(`I69-3` =="1" & Norway =="0" & OSKA2016 == "0") %>% add_count()
Osa_only <- match_all %>% filter(`I69-3` =="0" & Norway =="0" & OSKA2016 == "1") %>% add_count()
Nor_only <- match_all %>% filter(`I69-3` =="0" & Norway =="1" & OSKA2016 == "0") %>% add_count()



# charles
Oki <- venn_diagram$`I69-3`
Oki <- Oki[!is.na(Oki)]  # Remove NAs the old way
Nor <- venn_diagram$Norway  # This one has no NAs
Osa <- venn_diagram$OSKA2016 %>% discard(is.na)  # Remove NAs in fancier way.
l <- list(Oki=Oki, Osa=Osa, Nor=Nor)
sapply(l, length)
venn <- VennDiagram::get.venn.partitions(l, keep.elements = FALSE)
venn
plot(euler(l), shape = "ellipse", quantities = TRUE)
