################################################################################
### script to retrieve mutant identities for all phenotypes strains          ###
################################################################################

# load reference sequences and their information
projects <- read_csv("./data/global_key.csv")
refs <- read_csv("./data/Genetic_data/rpoB_reference_sequence.csv")
seqs <- readDNAStringSet("./data/Genetic_data/rpoB_references.fasta")
coordinates <- ALJEbinf::getAllCoordinates(seqs, "rpoB_Escherichia_coli_MG1655")

# load and complete reported mutations from Google sheet:
gs4_deauth() #function to inactive Google authorization
muts <- read_sheet("https://docs.google.com/spreadsheets/d/1-OihBy1WlfU46_npueE-ZJhmo1eH4eJf7BslsAVh6jc/edit?usp=sharing",  col_types = "c") #load google sheet
muts <- muts |>
  filter(substr(Ref_code, 1, 7) == "ALJElab") |>
  filter(Species %in% projects$Species) |>
  ALJEbinf::fillMutationsTable(refs, seqs, coordinates) # generate a new data frame using Google sheet, reference sequences and coordinates   

write_csv(muts, "./data/Genetic_data/muts.csv")
