FROM rocker/rstudio:latest

# Give rstudio user passwordless sudo access
RUN echo "rstudio ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers.d/rstudio

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    libxml2-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libglpk-dev \
    libfreetype6-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libgit2-dev \
    libbz2-dev \
    liblzma-dev \
    zlib1g-dev \
    libsodium-dev \
    && rm -rf /var/lib/apt/lists/*

# Configure R library path
RUN echo '.libPaths(c("/usr/local/lib/R/site-library", .libPaths()))' >> /home/rstudio/.Rprofile

# Set up RStudio config directories
RUN mkdir -p /home/rstudio/.config/rstudio/themes

# Ensure rstudio owns their home directory files
RUN chown -R rstudio:rstudio /home/rstudio