#create directory to download packages
brainArray_dir <- "brainArray"

if (!dir.exists(brainArray_dir))
  dir.create(brainArray_dir)

# Brain array download sources
# http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/25.0.0/entrezg.asp
# http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/CDF_download.asp

#Download brainArray custom CDF packages for each platform
install_brain_array <- function(id) {
  download.file(paste0("http://mbni.org/customcdf/25.0.0/entrezg.download/", id, "_25.0.0.tar.gz"),
                paste0("brainArray/", id, "_25.0.0.tar.gz"))
  install.packages(paste0("brainArray/", id, "_25.0.0.tar.gz"), repos = NULL, type = "source")
  unlink(paste0("brainArray/", id, "_25.0.0.tar.gz"))
}

#install brainArray custom CDF packages
install_brain_array("hugene10sthsentrezgprobe")   # Affymetrix Human Gene 1.0 ST Array [transcript (gene) version]
install_brain_array("huex10sthsentrezgprobe")     # Affymetrix Human Exon 1.0 ST Array [transcript (gene) version]
install_brain_array("u133aaofav2hsentrezgprobe")  # Affymetrix GeneChip HT-HG_U133A Early Access Array
install_brain_array("hgu95av2hsentrezgprobe")     # Affymetrix Human Genome U95 Version 2 Array
install_brain_array("hgu133ahsentrezgprobe")      # Affymetrix Human Genome U133A Array
install_brain_array("hgu133a2hsentrezgprobe")     # Affymetrix Human Genome U133A 2.0 Array
install_brain_array("hgu133plus2hsentrezgprobe")  # Affymetrix Human Genome U133 Plus 2.0 Array
install_brain_array("hgu133bhsentrezgprobe")      # Affymetrix Human Genome U133B Array


#Install annotation packages
BiocManager::install(c("pd.hugene.1.0.st.v1", "pd.huex.1.0.st.v2", "pd.ht.hg.u133a", "pd.hg.u133a", "pd.hg.u133a.2", "pd.hg.u133.plus.2", "pd.hg.u133b", "pd.hg.u95av2"))

# load brainArray libraries
library("hugene10sthsentrezgprobe")
library("huex10sthsentrezgprobe")
library("u133aaofav2hsentrezgprobe")
library("hgu133ahsentrezgprobe")
library("hgu133a2hsentrezgprobe")
library("hgu133plus2hsentrezgprobe")
library("hgu133bhsentrezgprobe")
library("hgu95av2hsentrezgprobe")

#load annotation package libraries
library("pd.hugene.1.0.st.v1")
library("pd.huex.1.0.st.v2")
library("pd.ht.hg.u133a")
library("pd.hg.u133a")
library("pd.hg.u133a.2")
library("pd.hg.u133.plus.2")
library("pd.hg.u133b")
library("pd.hg.u95av2")

#pd.hg.u95av2	Platform Design Info for Affymetrix HG_U95Av2
#hgu95av2.db	Affymetrix Human Genome U95 Set annotation data (chip hgu95av2)
