FROM mcr.microsoft.com/devcontainers/cpp:1-ubuntu-20.04

ARG REINSTALL_CMAKE_VERSION_FROM_SOURCE="3.22.2"

# Optionally install the cmake for vcpkg
COPY ./reinstall-cmake.sh /tmp/

RUN if [ "${REINSTALL_CMAKE_VERSION_FROM_SOURCE}" != "none" ]; then \
        chmod +x /tmp/reinstall-cmake.sh && /tmp/reinstall-cmake.sh ${REINSTALL_CMAKE_VERSION_FROM_SOURCE}; \
    fi \
    && rm -f /tmp/reinstall-cmake.sh

# [Optional] Uncomment this section to install additional vcpkg ports.
# RUN su vscode -c "${VCPKG_ROOT}/vcpkg install <your-port-name-here>"

# Use apt mirror site
RUN sed -i -e "s/http:\/\/archive\.ubuntu\.com\/ubuntu\//http:\/\/ubuntu.cs.nctu.edu.tw\/ubuntu\//" /etc/apt/sources.list

# [Optional] Uncomment this section to install additional packages.
RUN apt-get update && export DEBIAN_FRONTEND=noninteractive \
    && apt-get -y install --no-install-recommends gcc-10 clang-11 \
    && apt-get -y install --no-install-recommends linux-tools-$(uname -r) \
    && apt-get -y install --no-install-recommends tmux \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*