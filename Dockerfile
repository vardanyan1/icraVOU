# Use a lightweight base image
FROM ubuntu:latest

# Set the maintainer
LABEL maintainer="Vazgen Vardanyan <your.email@example.com>"

# Set the environment variables
ENV BASHRC=/root/.bashrc \
    PGPLOT_DIR=/pgplot \
    PGPLOT_DEV=/cps \
    RA=0 \
    DEC=0 \
    UUID=""

# Update and install dependencies
RUN apt-get update && apt-get install -y \
    git \
    gcc \
    gfortran \
    wget \
    curl \
    libpng-dev \
    python2.7 \
    python2-dev \
    python3 \
    python3.10 \
    unzip \
    libffi-dev \
    libssl-dev \
    python3-pip

# Manually install an older version of pip that supports Python 2
RUN curl https://bootstrap.pypa.io/pip/2.7/get-pip.py --output get-pip.py \
    && python2.7 get-pip.py

# Install PyYAML that is compatible with Python 2.7
RUN pip2 install 'pyyaml<5.4'

# Install EADA
RUN pip2 install https://github.com/chbrandt/eada/archive/0.9.7.5.zip

RUN pip3 install requests bs4 timeout_decorator lxml

RUN pip3 install  pandas astropy astroquery

# Create necessary directories for EADA
RUN mkdir -p /root/.config/eada

# Install PGPlot
RUN apt-get install -y pgplot5 libpq-dev

# Install dependencies
RUN pip2 install psycopg2-binary hsluv==5.0.0 pandas

# Clone VOU_Blazars repository
RUN git clone https://github.com/ecylchang/VOU_Blazars.git

# Copy all fortran files to the specific directory
COPY ./vou_new_files/fort_files/ /VOU_Blazars/bin/fort/

# Copy all other files to another directory
COPY ./vou_new_files/bin_files/ /VOU_Blazars/bin/

# Set the working directory
WORKDIR /VOU_Blazars/bin/fort

# Modify compile.sh to remove -mcmodel=medium flag
RUN sed -i 's/-mcmodel=medium//' compile.sh

# Compile mylib.f
RUN gfortran -c mylib.f -ffixed-line-length-500

# Compile other Fortran programs
RUN ./compile.sh

# Set the working directory
WORKDIR /VOU_Blazars

# Change the permission of the script to be executable
RUN chmod +x ./bin/vou-blazars-hybrid.sh

# Copy the Python script to the container
COPY app.py .
COPY ztf .
COPY conesearch_files ./conesearch_files

# Define a volume for the /work_dir directory
VOLUME /work_dir

# Set the entrypoint command to run the Python script
ENTRYPOINT ["python2", "app.py"]
