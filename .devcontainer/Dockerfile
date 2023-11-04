FROM mcr.microsoft.com/devcontainers/cpp:1-ubuntu-20.04

ARG PERSISTENCE_DIR="/persistence"
ARG BASH_HISTORY="${PERSISTENCE_DIR}/.bash_history"

WORKDIR /tmp

# Remove default vscode user and replace apt mirror site
RUN deluser vscode --remove-all-files \
    && sed -i -e "s/http:\/\/archive\.ubuntu\.com\/ubuntu\//http:\/\/ubuntu.cs.nctu.edu.tw\/ubuntu\//" /etc/apt/sources.list

RUN apt-get update && export DEBIAN_FRONTEND=noninteractive \
    && apt-get -y install --no-install-recommends gcc-10 g++-10 clang-11 \
    && apt-get -y install --no-install-recommends linux-tools-$(uname -r) \
    && apt-get -y install --no-install-recommends tmux imagemagick \
    && apt-get autoremove -y \
    && apt-get clean -y \
    && rm -rf /var/lib/apt/lists/*

# Set gcc-10, g++-10, clang-11 as default
RUN rm /usr/bin/gcc /usr/bin/g++ /usr/bin/gcov /usr/bin/clang /usr/bin/clang++ \
    && ln -s /usr/bin/gcc-10 /usr/bin/gcc \
    && ln -s /usr/bin/g++-10 /usr/bin/g++ \
    && ln -s /usr/bin/gcov-10 /usr/bin/gcov \
    && ln -s /usr/bin/clang-11 /usr/bin/clang \
    && ln -s /usr/bin/clang++-11 /usr/bin/clang++

# Install hyperfine
RUN wget --progress=dot:giga https://github.com/sharkdp/hyperfine/releases/download/v1.18.0/hyperfine-musl_1.18.0_amd64.deb \
    && dpkg -i hyperfine-musl_1.18.0_amd64.deb \
    && rm hyperfine-musl_1.18.0_amd64.deb

# Install TerminalImageViewer (tiv)
RUN git clone --depth 1 -b v.1.2 https://github.com/stefanhaustein/TerminalImageViewer.git
WORKDIR /tmp/TerminalImageViewer/src
RUN make install \
    && rm -rf /tmp/TerminalImageViewer
WORKDIR /tmp

# Persist bash history
RUN SNIPPET="export PROMPT_COMMAND='history -a' && export HISTFILE=${BASH_HISTORY}" \
    && echo "$SNIPPET" >> "$HOME/.bashrc"