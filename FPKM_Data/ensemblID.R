library(biomaRt)
#Not all ids resulted in a gene name


#get ensembl ids from expression matrix
ensembl.ids <- expressionMatrix[,1]

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

#replaces gene ids with names
expressionMatrix <- expressionMatrix[order(expressionMatrix$gene_id),]

expressionMatrix[,1] <- gene_names[,2]
