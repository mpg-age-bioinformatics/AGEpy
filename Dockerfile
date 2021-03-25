FROM python:3.8-slim

RUN apt-get update && apt-get install -yq --no-install-recommends git gcc g++ libz-dev imagemagick imagemagick-doc && apt-get clean && rm -rf /var/lib/apt/lists/*

RUN pip3 install git+https://github.com/mpg-age-bioinformatics/AGEpy.git
