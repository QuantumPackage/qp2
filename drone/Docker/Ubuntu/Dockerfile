ARG UBUNTU_VERSION=20.04
FROM ubuntu:${UBUNTU_VERSION} AS builder

# Timezone for tzdata
ARG tz=Etc/UTC
RUN echo $tz > /etc/timezone && rm -rf /etc/localtime

# Install
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y \
 git \
 curl \
 wget \
 python3 \
 gfortran \
 gcc \
 g++ \
 make \
 build-essential \
 rsync \
 unzip \
 libopenblas-dev \
 pkg-config \
 m4

RUN ln -s /usr/bin/python3 /usr/bin/python


