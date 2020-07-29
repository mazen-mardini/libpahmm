# ==============================================================================
#  libpahmm - library for paHMM-Tree, a phylogenetic tree estimator
#
#  Copyright (c) 2020 Mazen Mardini.
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#      This program is distributed in the hope that it will be useful,
#      but WITHOUT ANY WARRANTY; without even the implied warranty of
#      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#      GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses>.
# ==============================================================================

from cffi import FFI
import sys
import pcpp
from io import StringIO


def get_preprocessed_headers() -> str:
    p = pcpp.Preprocessor()
    p.add_path("auxiliary/fake_c_headers")
    p.define("PAHMM_NO_EXPORT_IMPORT")
    p.line_directive = None
    with open("include/cpahmm.h") as header_file:
        p.parse(header_file.read())
    output = StringIO()
    p.write(output)
    return output.getvalue()


ffibuilder = FFI()
ffibuilder.cdef(get_preprocessed_headers())

# set_source() gives the name of the python extension module to
# produce, and some C source code as a string.  This C code needs
# to make the declarated functions, types and globals available,
# so it is often just the "#include".
if sys.platform.startswith('linux'):
    linker_args = "-Wl,-rpath=$ORIGIN/pahmm,--enable-new-dtags"
elif sys.platform.startswith('darwin'):
    linker_args = "-Wl,-rpath,@loader_path/pahmm"
else:
    linker_args = ""

ffibuilder.set_source("_pahmm_cffi",
                      """
                           #include "cpahmm.h"   // the C header of the library
                      """,
                      extra_link_args=[linker_args],
                      library_dirs=["./python/pahmm"],
                      include_dirs=["./include"],
                      libraries=["pahmm"])

if __name__ == "__main__":
    ffibuilder.compile(verbose=True, tmpdir="./build/python")
