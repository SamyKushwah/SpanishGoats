#path to where you want data to be stored
data_path <- "S:/CGS/data"

# Create the data folder if it doesn't exist
if (!dir.exists(data_path)) {
  dir.create(data_path)
}


#This is the process to be repeated for each "data" table
#Need to figure out how to reinsert them back into the table, order is not maintained
#Not all ids resulted in a gene name

#path to text file with ensembl ids
ids_path <- "S:/CGS/data/Data1.txt"

#store ensembl ids in ensembl.eds **must be a .txt file with a header (column name)
ensembl.ids <- read.delim(ids_path, header = T)

#lists ensembl databases, we use genes
listEnsembl()

ensembl <- useEnsembl(biomart = "genes")

#list of ensembl data sets, filtered for "goat"
datasets <- listDatasets(ensembl)

#connection with the database
ensembl.con <- useMart("ensembl", dataset = 'chircus_gene_ensembl', )

#list of attributes (we want gene name)
attr <- listAttributes(ensembl.con)

#list of filters (we are using the "gene stable")
filters <- listFilters(ensembl.con)

#saves the results of the query into gene_names
gene_names <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                    filters = "ensembl_gene_id",
                    values = ensembl.ids$gene_id,
                    mart = ensembl.con)
