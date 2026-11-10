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

# Rust toolchain for the tool layer (`rust/`, pinned in
# rust/rust-toolchain.toml). `./configure` only checks that `cargo` exists;
# the build itself runs `cargo build --release` through
# scripts/compilation/build_rust_tools.sh (the `rust_target` rule).
ENV RUSTUP_HOME=/usr/local/rustup \
    CARGO_HOME=/usr/local/cargo \
    PATH=/usr/local/cargo/bin:$PATH
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs \
      | sh -s -- -y --profile minimal --default-toolchain none && \
    rustup toolchain install 1.92.0 --profile minimal


