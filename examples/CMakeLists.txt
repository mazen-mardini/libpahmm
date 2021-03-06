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

project(examples)

add_executable(examples
    main.cpp
)

target_link_libraries(examples PUBLIC pahmm)
target_link_libraries(examples PUBLIC paHMM-dist)
target_link_libraries(examples PUBLIC pthread)
target_include_directories(examples PRIVATE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../src>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../dlib>
)
