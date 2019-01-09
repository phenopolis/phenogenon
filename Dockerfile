FROM ubuntu:latest

LABEL maintainer="Jing Yu <logust79@gmail.com>"
LABEL Description="phenogenon" Version="1.0.0"

ENV USER root
ENV TMPDIR /tmp

ENV TZ=GMT
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt-get update
RUN apt-get install -y --no-install-recommends build-essential checkinstall \
  libreadline-gplv2-dev libncursesw5-dev libssl-dev libsqlite3-dev tk-dev \
  libgdbm-dev libc6-dev libbz2-dev zlib1g-dev openssl libffi-dev python \
  python-dev python-setuptools python-pip python-numpy wget git tabix jq bcftools
RUN rm -rf /var/lib/apt/lists/*

# Add source code
ADD . /phenogenon
WORKDIR /phenogenon

RUN python setup.py install
RUN pip install --no-cache-dir coveralls

CMD phenogenon -h
