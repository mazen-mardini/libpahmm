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

__ALL__ = ["PAHMMError", "BandingEstimator", "Sequences"]

from _pahmm_cffi import lib as _lib, ffi as _ffi
from typing import AnyStr, Union
from pathlib import Path


class PAHMMError(Exception):
    """PaHMM exception class.

    Errors coming from pahmm will generate an exception of this type.
    """

    def __init__(self, message, cffi_be):
        c_error = _ffi.string(_lib.ebc_be_last_error_msg(cffi_be)).decode("utf8")
        super(PAHMMError, self).__init__(message + " " + c_error if c_error else "")


class BandingEstimator:
    """Loads and generates Sequences-objects for MSA (distance-) calculations.

    Use this class to load a file or string, then use any of the desired execute_*-method
    to generate a Sequences-object. The Sequences-object can then be used to access the
    loaded sequences, their names and distances between any given sequences.
    """

    def __init__(self):
        self.__be = _lib.ebc_be_create()

        # Tell the garbage collector how to free the resources
        self.__be = _ffi.gc(self.__be, _lib.ebc_be_free)

    @staticmethod
    def _path_to_bytes(path: Union[Path, AnyStr]) -> bytes:
        """Used to convert strings and bytearrays to bytes-objects.
        """
        if isinstance(path, bytearray):
            return bytes(path)
        elif isinstance(path, str):
            return path.encode("utf8")
        elif isinstance(path, Path):
            return str(path).encode("utf8")

    def execute_gtr_model(self, arg1: Union[None, float] = None, arg2: Union[None, float] = None,
                          arg3: Union[None, float] = None, arg4: Union[None, float] = None,
                          arg5: Union[None, float] = None) -> Union[None, "Sequences"]:
        """Generate a Sequences-object that will use the GTR model for distance calculations.

        This model can only be used for nucleotides.
        """

        if None in (arg1, arg2, arg3, arg4, arg5):
            # One if the arguments is None, choose to estimate the model parameters instead.
            sequences = _lib.ebc_be_execute_gtr_modelv2(self.__be)
        else:
            # Use the supplied model parameters
            sequences = _lib.ebc_be_execute_gtr_model(self.__be, arg1, arg2, arg3, arg4, arg5)

        if sequences == _ffi.NULL:
            raise PAHMMError("Could not execute GTR model.", self.__be)

        # Tell the garbage collector how to free the resources
        sequences = _ffi.gc(sequences, _lib.ebc_seq_free)

        return Sequences(sequences)

    def execute_hky85_model(self, arg: Union[None, float] = None) -> Union[None, "Sequences"]:
        """Generate a Sequences-object that will use the HKY85 model for distance calculations.

        This model can only be used for nucleotides.
        """

        if arg is None:
            # The argument is None, choose to estimate the model parameters instead.
            sequences = _lib.ebc_be_execute_hky85_modelv2(self.__be)
        else:
            # Use the supplied model parameters
            sequences = _lib.ebc_be_execute_hky85_model(self.__be, arg)

        if sequences == _ffi.NULL:
            raise PAHMMError("Could not execute HKY85 model.", self.__be)

        # Tell the garbage collector how to free the resources
        sequences = _ffi.gc(sequences, _lib.ebc_seq_free)

        return Sequences(sequences)

    def execute_jtt_model(self) -> Union[None, "Sequences"]:
        """Generate a Sequences-object that will use the JTT model for distance calculations.

        This model can only be used for amino-acids.
        """

        sequences = _lib.ebc_be_execute_jtt_model(self.__be)

        if sequences == _ffi.NULL:
            raise PAHMMError("Could not execute JTT model.", self.__be)

        # Tell the garbage collector how to free the resources
        sequences = _ffi.gc(sequences, _lib.ebc_seq_free)

        return Sequences(sequences)

    def execute_lg_model(self) -> Union[None, "Sequences"]:
        """Generate a Sequences-object that will use the LG model for distance calculations.

        This model can only be used for amino-acids.
        """

        sequences = _lib.ebc_be_execute_lg_model(self.__be)

        if sequences == _ffi.NULL:
            raise PAHMMError("Could not execute LG model.", self.__be)

        # Tell the garbage collector how to free the resources
        sequences = _ffi.gc(sequences, _lib.ebc_seq_free)

        return Sequences(sequences)

    def execute_wag_model(self) -> Union[None, "Sequences"]:
        """Generate a Sequences-object that will use the WAG model for distance calculations.

        This model can only be used for amino-acids.
        """

        sequences = _lib.ebc_be_execute_wag_model(self.__be)

        if sequences == _ffi.NULL:
            raise PAHMMError("Could not execute WAG model.", self.__be)

        # Tell the garbage collector how to free the resources
        sequences = _ffi.gc(sequences, _lib.ebc_seq_free)

        return Sequences(sequences)

    def set_indel_paramters(self, nb_probability: Union[None, float] = None,
                            rate: Union[None, float] = None):
        """Set the nb-probability and rate parameters.
        """

        if None in (nb_probability, rate):
            # One of the arguments was None, choose to estimate the parameters instead.
            _lib.ebc_be_unset_indel_parameters(self.__be)
        else:
            _lib.ebc_be_set_indel_parameters(self.__be, nb_probability, rate)

    def __getattr__(self, key):
        """Get general attributes for this banding estimator.

        There are two general attributes: 'alpha' and 'gamma_rate_categories'
        """

        if key == "alpha":
            return self.__be[0].alpha if not self.__be[0].estimate_alpha else None
        elif key == "gamma_rate_categories":
            return self.__be[0].gamma_rate_categories if not self.__be[0].estimate_categories else None

    def __setattr__(self, key, value):
        if key == "alpha":
            if value is None:
                _lib.ebc_be_unset_alpha(self.__be)
            else:
                _lib.ebc_be_set_alpha(self.__be, value)
        elif key == "gamma_rate_categories":
            if value is None:
                _lib.ebc_be_unset_categories(self.__be)
            else:
                _lib.ebc_be_set_categories(self.__be, value)
        else:
            super().__setattr__(key, value)

    def set_str_input(self, fasta: Union[str, bytes, bytearray]):
        """Decide to load input from a string/bytes/byterarray.

        For better performance use set_file_input instead.
        """

        result: bool = \
            _lib.ebc_be_set_input(self.__be, fasta.encode("utf8") if isinstance(fasta, str) else fasta)

        if not result:
            raise PAHMMError("Could not read FASTA input.", self.__be)

    def set_file_input(self, fasta_path: Union[Path, AnyStr]):
        """Decide to load input from a file.

        This is faster than set_str_input because all data will be read by libpahmm directly,
        bypassing Python.
        """

        fasta_path = self._path_to_bytes(fasta_path)
        result: bool = _lib.ebc_be_set_input_from_file(self.__be, fasta_path)

        if not result:
            raise PAHMMError("Could not read FASTA file input.", self.__be)


class Sequences:
    """A class used to access sequences and calculate distances between them.

    Note: Do not instantiate an object of this class directly. Use one of
    BandingEstimator's execute_*-methods for that.
    """

    def __init__(self, seq):
        self._seq_count = _lib.ebc_seq_count(seq)
        self._seq = seq

    def __len__(self):
        """Get the number of sequenes.
        """
        return self._seq_count

    def __getitem__(self, seq_id):
        """Get a sequence using a sequence number/ID.
        """
        return _ffi.string(_lib.ebc_seq_get_sequence(self._seq, seq_id))

    def get_distance(self, seq_id1: int, seq_id2: int):
        """Retrieve a distance between two sequences using their numbers/IDs.
        """
        return _lib.ebc_seq_get_distance(self._seq, seq_id1, seq_id2)

    def get_distance_from_names(self, seq_name1: Union[bytes, bytearray], seq_name2: Union[bytes, bytearray]):
        """Retrieve a distance between two sequences using their names.
        """
        return _lib.ebc_seq_get_distance_from_names(self._seq, seq_name1, seq_name2)

    def get_seq_name(self, seq_id: int):
        """Get a sequence using its number or ID.
        """
        return _ffi.string(_lib.ebc_seq_get_name(self._seq, seq_id))

    def get_sequence_from_name(self, seq_name: Union[bytes, bytearray]):
        """Get a sequence using its name.
        """
        return _ffi.string(_lib.ebc_seq_get_sequence_from_name(self._seq, seq_name))

