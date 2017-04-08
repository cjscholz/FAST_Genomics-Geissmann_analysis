############################################################ 
# Dockerfile to run the Seurat preprocessing and clustering
# of single-cell RNA-seq data in FAST Genomics format
############################################################ 

# File Author / Maintainer 
MAINTAINER Claus J. Scholz <cscholz@uni-bonn.de>


FROM bioconductor/release_core2
ADD ./install_scripts /opt/install_scripts
RUN Rscript /opt/install_scripts/install_Seurat.R
RUN apt-get -y install default-jre

ADD ./scripts /opt/scripts
CMD ["Rscript", "--default-packages=base,utils,grDevices,graphics,stats,methods", "/opt/scripts/raw2cluster.R"]
