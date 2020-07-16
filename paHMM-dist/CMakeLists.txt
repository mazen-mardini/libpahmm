##==============================================================================
## libpahmm - library for paHMM-Tree, a phylogenetic tree estimator
##
## Copyright (c) 2020 Mazen Mardini.
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
##     This program is distributed in the hope that it will be useful,
##     but WITHOUT ANY WARRANTY; without even the implied warranty of
##     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##     GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses>.
##==============================================================================

project(paHMM-dist LANGUAGES CXX)

cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

file(GLOB SOURCE_PATHS
    "src/core/*.h"
    "src/core/*.cpp"
    "src/heuristics/*.h"
    "src/heuristics/*.cpp"
    "src/hmm/*.h"
    "src/hmm/*.cpp"
    "src/models/*.h"
    "src/models/*.cpp"
)

# Settings for dlib
# These settings will disable unused features. You may reintroduce these features,
# but remember that this might imply more external dependencies which might reduce
# portability.
set(DLIB_PNG_SUPPORT OFF)
set(DLIB_GIF_SUPPORT OFF)
set(DLIB_JPEG_SUPPORT OFF)
set(DLIB_USE_BLAS OFF)       # Might impact performance, to use BLAS reenable this.
set(DLIB_USE_LAPACK OFF)     # Might impact performance, to use LAPACK reenable this.
set(DLIB_USE_MKL_FFT OFF)    # Might impact performance, to use MKL_FFT reenable this.
set(DLIB_USE_CUDA OFF)       # Might impact performance, to use CUDA reenable this (Nvidia GPU only).
set(DLIB_NO_GUI_SUPPORT ON)
# This MUST be on:
set(DLIB_DISABLE_ASSERTS ON)

# Include dlib
add_subdirectory(dlib)

# Settings for dlib
set_property(TARGET dlib PROPERTY CXX_STANDARD 11)
#add_definitions(-DDLIB_DISABLE_ASSERTS)

add_library(paHMM-dist STATIC ${SOURCE_PATHS} dlib/dlib/all/source.cpp)
target_link_libraries(paHMM-dist PRIVATE dlib)



target_include_directories(paHMM-dist PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>/src
)

set_property(TARGET dlib PROPERTY POSITION_INDEPENDENT_CODE ON)
#set_property(TARGET paHMM-dist PROPERTY POSITION_INDEPENDENT_CODE ON)
