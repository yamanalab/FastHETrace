#!/bin/bash

set -eu

cmake -S . -B build-normal \
        -DCMAKE_EXPORT_COMPILE_COMMANDS=yes \
        -DCMAKE_INSTALL_PREFIX=$HOME/.local \
        -DWITH_INTEL_HEXL=no
cmake -S . -B build-hexl \
        -DCMAKE_EXPORT_COMPILE_COMMANDS=yes \
        -DCMAKE_INSTALL_PREFIX=$HOME/.local \
        -DWITH_INTEL_HEXL=yes
cmake --build build-normal -j
cmake --build build-hexl -j
