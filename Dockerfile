FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    git \
    pkg-config \
    python3-dev \
    python3-pip \
    python3-numpy \
    && rm -rf /var/lib/apt/lists/*

RUN git clone --recurse-submodules https://github.com/xlnt-community/xlnt.git /tmp/xlnt \
    && cd /tmp/xlnt \
    && mkdir build && cd build \
    && cmake -D STATIC=ON .. \
    && make -j$(nproc) \
    && make install \
    && rm -rf /tmp/xlnt

WORKDIR /workspace

COPY . /workspace

RUN make clean && make


# Default command to run the stroke analysis program with specified input and output files
CMD ["./stroke_analysis", "--input", "data/Canoe_xsens_dot.csv", "--output", "results"]
