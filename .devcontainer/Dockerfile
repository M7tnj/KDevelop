FROM rocker/rstudio:latest


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


RUN echo '.libPaths(c("/usr/local/lib/R/site-library", .libPaths()))' >> /home/rstudio/.Rprofile


RUN mkdir -p /home/rstudio/.config/rstudio/themes


RUN chown -R rstudio:rstudio /home/rstudio