FROM ubuntu:16.04
RUN apt-get update && apt-get install -y python sudo dstat
RUN apt-get install -y python-pip wget zlib1g-dev
RUN pip install -U pip
RUN pip install -U setuptools
RUN pip install numpy scipy pysam cython
RUN pip install h5py pandas ArgumentParser
RUN pip install lxml

#copy contents of tasks/<taskname>/build/src on the build host into /opt/src on the docker
COPY src/ /opt/src/

WORKDIR /opt/src
