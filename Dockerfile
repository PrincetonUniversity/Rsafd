FROM r-base:latest

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libpng-dev \
    libjpeg-dev \
    gfortran \
    && rm -rf /var/lib/apt/lists/*

# Copy the R package source
COPY . /usr/src/Rsafd

# Set working directory
WORKDIR /usr/src/Rsafd

# Install dependencies
RUN Rscript -e "install.packages(c('timeDate', 'quadprog', 'quantreg', 'plot3D', 'robustbase', 'scatterplot3d', 'splines', 'tseries', 'glasso', 'qgraph', 'reticulate', 'keras', 'rgl', 'glmnet'), repos='https://cran.rstudio.com')"

# Install package without test load (skips namespace load check)
RUN R CMD INSTALL --no-test-load .

# Explicit test load (will fail here if namespace truly broken)
RUN Rscript -e "library(Rsafd)"
