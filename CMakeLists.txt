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

cmake_minimum_required(VERSION 3.10)

# Set project name and version
project(libpahmm LANGUAGES CXX VERSION 0.1.0)

# Set minimum C++ standard to C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Define options
option(BUILD_STATIC "Build libpahmm as static" OFF)

# The paHMM library
list(APPEND SOURCE_FILES
    include/cpahmm.h
    src/cpahmm.cpp
    src/cpahmm_p.h
    src/cpahmm_p.cpp
    src/StreamParser.hpp
    src/StreamParser.cpp
    README.md
    LICENSE

    dlib/dlib/all/source.cpp
)

if (BUILD_STATIC)
    add_library(pahmm STATIC ${SOURCE_FILES})
    set_property(TARGET pahmm PROPERTY POSITION_INDEPENDENT_CODE ON)
else ()
    add_library(pahmm SHARED ${SOURCE_FILES})
endif ()

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

# Add dependencies
add_subdirectory(dlib)
add_subdirectory(paHMM-dist)

# Settings for dlib
set_property(TARGET dlib PROPERTY CXX_STANDARD 11)
#add_definitions(-DDLIB_NO_GUI_SUPPORT)
# Because something is turning it on, we have to turn it off:
#add_definitions(-DDLIB_DISABLE_ASSERTS)

target_link_libraries(pahmm PRIVATE paHMM-dist dlib)

target_include_directories(pahmm PRIVATE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)
target_include_directories(pahmm PUBLIC
    "${PROJECT_BINARY_DIR}"
)

target_compile_definitions(pahmm PRIVATE PAHMM_LIBRARY)

# Installation
install(TARGETS pahmm DESTINATION lib)
install(FILES include/cpahmm.h DESTINATION include)

if (UNIX)
    add_custom_target(libpahmm.pc ALL
            ${CMAKE_COMMAND} -DIN="${CMAKE_CURRENT_SOURCE_DIR}/libpahmm.pc.in" -DOUT="libpahmm.pc"
            -DPREFIX="${CMAKE_INSTALL_PREFIX}" -DVERSION="${libpahmm_VERSION}"
            -P "${CMAKE_CURRENT_SOURCE_DIR}/pkgconfig.cmake"
            COMMENT "Creating pkg-config file")

    install(FILES "${CMAKE_CURRENT_BINARY_DIR}/libpahmm.pc"
            DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/pkgconfig")
endif ()

# Setup the examples
add_subdirectory(examples)


