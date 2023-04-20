# Download base image ubuntu 16.04
# We choose this old image to install an old version of Maucalay2 (BNReduction script)
FROM ubuntu:16.04

# LABEL about the custom image
LABEL maintainer="bblodfon@gmail.com"
LABEL version="1.0"
LABEL description="This is custom Docker Image for the druglogics-synergy module."

# Install Java 8
RUN apt-get update && \
	apt-get install -y openjdk-8-jdk && \
	apt-get clean && \
	rm -rf /var/lib/apt/lists/* && \
	rm -rf /var/cache/oracle-jdk8-installer;

# Fix certificate issues, found as of 
# https://bugs.launchpad.net/ubuntu/+source/ca-certificates-java/+bug/983302
RUN apt-get update && \
	apt-get install -y ca-certificates-java && \
	apt-get clean && \
	update-ca-certificates -f && \
	rm -rf /var/lib/apt/lists/* && \
	rm -rf /var/cache/oracle-jdk8-installer;

# Setup JAVA_HOME, this is useful for docker commandline                        
ENV JAVA_HOME /usr/lib/jvm/java-8-openjdk-amd64/                                
RUN export JAVA_HOME

# get git and other commands needed to install BNReduction.sh
RUN apt-get update && \
	apt-get install -y git && \
	apt-get install -y wget && \
	apt-get install -y unzip && \
	apt-get install -y install-info && \
	apt-get install -y bzip2 && \
	apt-get install -y g++ && \
	apt-get install -y make;

# Install BNReduction.sh script
RUN git clone https://github.com/druglogics/druglogics-dep.git

RUN cd druglogics-dep && \
	dpkg -i dep/libpari-gmp3_2.5.0-2ubuntu1_amd64.deb && \
	dpkg -i dep/Macaulay2-1.6-common.deb && \
	# liblapack3 install problem
	apt-get install liblapack3 -y && \
	#apt-get -f install -y && \
	dpkg -i dep/Macaulay2-1.6-amd64-Linux-Ubuntu-14.04.deb && \
	unzip dep/bnet_reduction-master.zip -d dep && \
	tar jxfv dep/boost_1_55_0.tar.bz2 -C dep/bnet_reduction-master && \
	cd dep/bnet_reduction-master/ && \
	make clean && \
	make install && \
	./Testing_BNReduction.sh;

# BNReduction druglogics-specific installation
ENV BNET_HOME /druglogics-dep/dep/bnet_reduction-master
RUN cd druglogics-dep && \
	cp BNReduction.sh $BNET_HOME && \
	rm dep/boost*;

# Install python9 and mpbn
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-py39_4.10.3-Linux-x86_64.sh && \
	bash Miniconda3-py39_4.10.3-Linux-x86_64.sh -b && \
	ln -s /root/miniconda3/bin/python3 /usr/bin/python && \
	ln -s /root/miniconda3/bin/conda /usr/bin/conda && \
	# should be Python 3.9
	python --version && \
	# should be conda 4.10
	conda --version && \
	rm Miniconda3-py39_4.10.3-Linux-x86_64.sh && \
	conda install -y -c colomoto -c potassco mpbn=1.6 && \
	cd druglogics-dep && \
	python mpbn-attractors.py model.bnet;

ENV MPBN_HOME /druglogics-dep

# Get druglogics-synergy via git
#RUN git clone https://github.com/druglogics/druglogics-synergy.git && \
#	cd druglogics-synergy && \
#	git checkout v1.2.1

# copy Maven-build package .jar file
COPY target/synergy-1.2.1-jar-with-dependencies.jar /synergy-1.2.1-jar-with-dependencies.jar

# specify default command
ENTRYPOINT ["/usr/bin/java", "-cp", "/synergy-1.2.1-jar-with-dependencies.jar", "eu.druglogics.synergy.Launcher"]

