# CMake
cmake_minimum_required (VERSION 2.8)
# project
project (SAW)
# find all source code
# save to variable  DIR_SRCS
aux_source_directory(. DIR_SRCS)
# build target
add_executable(SAW ${DIR_SRCS})
