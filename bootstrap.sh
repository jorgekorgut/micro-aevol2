#!/bin/sh

set -eu

cd "$(dirname "$(realpath "$0")")"/boost

git submodule sync
git submodule update --init

wget 'https://boostorg.jfrog.io/artifactory/main/release/1.83.0/source/boost_1_83_0.tar.gz'
tar xvf boost_1_83_0.tar.gz
rm boost_1_83_0.tar.gz

ln -srf dynamic_bitset/include/boost/dynamic_bitset/dynamic_bitset.hpp \
        boost_1_83_0/boost/dynamic_bitset/dynamic_bitset.hpp
