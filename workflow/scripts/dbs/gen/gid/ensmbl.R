library(biomaRt)

# Parse args
orgms <- commandArgs(trailingOnly = TRUE)

get_gene_table <- function(dataset){
    
    # Connect to the Ensembl database
    ensembl <- useMart(
        biomart = 'ENSEMBL_MART_ENSEMBL',
        dataset = dataset,
        host    = 'https://www.ensembl.org'
    )

    # Specify the attributes to retrieve
    attributes <- c("ensembl_gene_id", "external_gene_name")
    # Retrieve the data
    gene_data <- getBM(
        attributes = attributes,
        mart = ensembl,
        useCache=FALSE,
        verbose=FALSE
    )
    colnames(gene_data) <- c('id', 'symbol')
    return(gene_data)
}

org_table <- list(
    'hg38'='hsapiens_gene_ensembl',
    'mm10'='mmusculus_gene_ensembl'
)

for (path_org in orgms) {
    org <- sub('^dbs/([^/]+)/.*$', '\\1', path_org)
    org <- org_table[org]
    gid <- get_gene_table(org)
    gid <- gid[gid$symbol != "", ]
    write.csv(x = gid, file = path_org, row.names=FALSE, quote=FALSE)
}
