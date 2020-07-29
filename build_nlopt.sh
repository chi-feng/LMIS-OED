#!/bin/bash

cd extern/nlopt
mkdir -p build
cd build
cmake ..
make -j8
sudo make install