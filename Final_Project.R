library(tidyverse)
library(ShortRead)
library(msa)
library(phangorn)
library(seqinr)
library(adegenet)
library(ape)
library(ggtree)
library(DECIPHER)
library(viridis)
library(ggplot2)
library(muscle)
library(purrr)


#   1. ALIGNING THE RAW SEQUENCES
#reading the fasta file and assigning it to a dataframe
raw_sequences <- readDNAStringSet("./ungulates.fasta")

#glimpse at the raw nucleotide sequences
raw_sequences

#Re-orienting nucleotide sequences 
oriented_sequences <- OrientNucleotides(raw_sequences)


#aligning the similarities between the sequences 
#and turning them all to the same length

#muscle alignment
muscle_aligned_sequences <- muscle(oriented_sequences)

muscle_aligned_sequences #glimpse at the muscle alignment

muscle_dnastr = as(muscle_aligned_sequences, "DNAStringSet") #this excises the non coding regions, and turns the multiple alignment into DNAstring

# change labels to show the Species names.
names(muscle_dnastr)

#We go into the location of the names, and then divide each string based on spaces.
#Then we select the second and third field. 
#We paste those fields together and separate them by a "_"
muscle_dnastr@ranges@NAMES <-
paste( muscle_dnastr@ranges@NAMES %>% str_split(pattern = " ") %>%map_chr(2),
       muscle_dnastr@ranges@NAMES %>% str_split(pattern = " ") %>%map_chr(3),
       sep = "_")

names(muscle_dnastr)

#viewing the sequences in a browser
BrowseSeqs(muscle_dnastr)



#We do this step in order to be able to build a tree
#writing the aligned data
writeXStringSet(muscle_dnastr, file ="ungulates_muscle_aligned.fasta")


#   2. BUILDING A TREE FROM THE ALIGNEMENT PRODUCED. 

#Reading the aligned file
dna <- read.alignment("ungulates_muscle_aligned.fasta", format = "fasta")

#Creating a Distance matrix for the aligned sequences 
d_matrix <- dist.alignment(dna, matrix = "similarity")
d_matrix

#Visualizing the distance matrix
#darker = farther, lighter= closer sequences
d_matrix_table <- as.data.frame(as.matrix(d_matrix))
table.paint(d_matrix_table, cleg=0, clabel.row = .5, clabel.col = .5)+
  scale_color_viridis()

tree <- nj(d_matrix)
class (tree)

#Ladderizing based on the distance matrix
tree <- ladderize(tree)

tree<- root(tree, outgroup = "Pan_troglodytes", resolve.root = TRUE)


names(muscle_dnastr)

plot(tree, cex = 0.6)

ggtree(tree, branch.length = "none") +
  #geom_tiplab()+
  geom_tiplab(hjust= -.3)+
  geom_text(aes(label=node), size=2.5, hjust=-.8)+
  geom_point2(aes(subset=(node==43)), shape=21, size=5, fill='red')+
  geom_point2(aes(subset=(node==59)), shape=21, size=5, fill='green')+
  geom_hilight(node=43, fill="orange") +
  geom_hilight(node=c(12,50), fill="blue") +
  geom_cladelabel(node=50, label="Zebras", 
                  color="purple", offset=5, align=TRUE) +
  geom_cladelabel(node=43, label="Donkeys", 
                  color="orange", offset=7, align=TRUE) +
  #geom_hilight(node=46, fill="gold")
  xlim(0,20)


