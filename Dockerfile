FROM ubuntu:latest

MAINTAINER Jing Yu <logust79@gmail.com>

LABEL Description="phenogenon" Version="1.0.0"

#ARG GITUSERNAME
#ARG GITPASSWORD

ENV USER root
ENV TMPDIR /tmp

RUN apt-get update -y

ENV TZ=GMT
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt-get install -y build-essential
RUN apt-get install -y checkinstall
RUN apt-get install -y libreadline-gplv2-dev
RUN apt-get install -y libncursesw5-dev
RUN apt-get install -y libssl-dev
RUN apt-get install -y libsqlite3-dev
RUN apt-get install -y tk-dev
RUN apt-get install -y libgdbm-dev
RUN apt-get install -y libc6-dev
RUN apt-get install -y libbz2-dev
RUN apt-get install -y zlib1g-dev
RUN apt-get install -y openssl
RUN apt-get install -y libffi-dev
RUN apt-get install -y python
RUN apt-get install -y python-dev
RUN apt-get install -y python-setuptools
RUN apt-get install -y python-pip
RUN apt-get install -y wget
RUN apt-get install -y git
RUN apt-get install -y tabix
RUN apt-get install -y bcftools

RUN pip install --upgrade pip
RUN pip install pymongo==3.4.0
RUN pip install biopython==1.68
RUN pip install pyliftover==0.3
RUN pip install bs4==0.0.1
RUN pip install scipy==0.19.0
RUN pip install pandas==0.20.2
RUN pip install pysam==0.11.2.2
RUN pip install fisher==0.1.5
RUN pip install coveralls




