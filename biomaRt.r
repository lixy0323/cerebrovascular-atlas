library(biomaRt)
    listMarts()
    umart = useMart('ensembl', host = "https://dec2021.archive.ensembl.org/")
    datalist  =listDatasets(umart)
    searchDatasets(mart = umart, pattern = "Human")
    searchDatasets(mart = umart, pattern = "Mouse")
    searchDatasets(mart = umart, pattern = "Zebrafish")
    mouse <- useMart('ensembl',dataset = "mmusculus_gene_ensembl",host = "https://dec2021.archive.ensembl.org/")
    zebrafish <- useMart('ensembl',dataset = "drerio_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
    human <- useMart('ensembl',dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
    if(F){
      g1 <- getBM(
        attributes = c("external_gene_name", "ensembl_gene_id","entrezid_gene_id"),
        mart = mouse,
        uniqueRows = T
      )
      
      gene.m2z <-
        getLDS(
          attributes = c("external_gene_name"),
          filters = "external_gene_name",
          values = genes_tobe_transfer(a_vecter),
          mart = mouse,
          attributesL = c("external_gene_name"),
          martL = zebrafish,
          uniqueRows = T
        )
    }
  }
