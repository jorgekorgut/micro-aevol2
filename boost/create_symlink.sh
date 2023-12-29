#!/bin/sh

# Get the path of the script
script_path=$(dirname $(realpath $0))

# Specify the target file and the symlink name
target=$script_path/boost_1_83_0/boost/dynamic_bitset/dynamic_bitset.hpp
symlink=$script_path/dynamic_bitset/include/boost/dynamic_bitset/dynamic_bitset.hpp

# Remove the target if it exists
if [ -L "$target" ]; then
    rm "$target"
fi

# Create the symlink
ln -s "$symlink" "$target" 