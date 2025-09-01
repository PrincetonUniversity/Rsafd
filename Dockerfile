FROM r-base:latest

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libpng-dev \
    libjpeg-dev \
    gfortran \
    && rm -rf /var/lib/apt/lists/*



# Install dependencies
RUN Rscript -e "install.packages(c('timeDate', 'quadprog', 'quantreg', 'plot3D', 'robustbase', 'scatterplot3d', 'splines', 'tseries', 'glasso', 'qgraph', 'reticulate', 'keras', 'rgl', 'glmnet'), repos='https://cran.rstudio.com')"

# Copy the R package source
COPY . /usr/local/lib/R/site-library/Rsafd

# Explicit test load (will fail here if namespace truly broken)
RUN R -q -e "library(Rsafd)"