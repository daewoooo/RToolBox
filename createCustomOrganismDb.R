## Creating custom human hg38 transcript database ##
####################################################
source("https://bioconductor.org/biocLite.R")
library(OrganismDbi)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

gd <- list(join1 = c(GO.db="GOID", org.Hs.eg.db="GO"), 
           join2 = c(org.Hs.eg.db = "ENTREZID", 
                     TxDb.Hsapiens.UCSC.hg38.knownGene = "GENEID"))

destination <- tempfile()
dir.create(destination)
makeOrganismPackage(pkgname = "Homo.sapiens.hg38", graphData = gd, 
                    organism = "Homo sapiens", version = "1.0.0", 
                    maintainer = "Maintainer<maintainer@email>", 
                    author = "DavidPorubsky", destDir = destination, 
                    license = "Artistic-2.0")

## >> Creating package in /tmp/RtmpmBnznk/file12a01043ad52/Homo.sapiens.hg38 

## installing package
install.packages('/tmp/RtmpmBnznk/file12a01043ad52/Homo.sapiens.hg38', repos = NULL, type="source")
                                                                                                                                                                                 organism = "Homo sapiens", version = "1.0.0", 
                                                                                                                                                                                 maintainer = "Maintainer<maintainer@email>", 
                                                                                                                                                                                 author = "Author Name", destDir = destination, 
                                                                                                                                                                                 license = "Artistic-2.0")#Creating package in C:/Users/arvis/AppData/Local/Temp/RtmpC8UxBO/file1fd4a713105/Homo.sapiens.hg38install.packages("C:/Users/arvis/AppData/Local/Temp/RtmpC8UxBO/file1fd4a713105/Homo.sapiens.hg38", repos = NULL, type="source")