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

# Install R dependencies and the package
RUN R -e "install.packages('remotes')"
RUN R -e "remotes::install_deps(dependencies = TRUE)"
RUN R -e "remotes::install_local(build = FALSE)"
