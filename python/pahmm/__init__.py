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

__all__ = ["PAHMMError", "BandingEstimator", "Sequences"]

from _pahmm_cffi import lib as _lib, ffi as _ffi
from typing import AnyStr, Union
from pathlib import Path


class PAHMMError(Exception):
    """PaHMM exception class.

    Errors coming from pahmm will generate an exception of this type.
    """

    def __init__(self, message, be=None):
        c_error = be.last_error_message() if be else None
        if c_error:
            super(PAHMMError, self).__init__(message + " " + c_error if c_error else "")
        else:
            super(PAHMMError, self).__init__(message)


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

    def apply_model(self, model):
        """
        After you have associated sequence data to the estimator (using
        set_str_input or set_file_input), use this function to apply a model to
        the data. Necessary parameters are estimated and rough distances are
        computed. This is a convenience function.
        """
        if model == 'HKY85':
            return self.execute_hky85_model()
        elif model == 'GTR':
            return self.execute_gtr_model()
        elif model == 'WAG':
            return self.execute_wag_model()
        elif model == 'LG':
            return self.execute_lg_model()
        elif model == 'JTT':
            return self.execute_jtt_model()
        else:
            raise Exception(f'Bug! The model "{model}" is not implemented by paHMM.')

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
        else:
            raise PAHMMError(f"Expected a path, got an object of type: {type(path).__name__}")

    def has_last_error(self):
        """Returns True if last error message is available, and False otherwise.
        """
        error_ptr = _lib.ebc_be_last_error_msg(self.__be)
        return error_ptr != _ffi.NULL

    def last_error_message(self) -> Union[None, str]:
        """Get last error message.

        Returns None if none is available.
        """
        error_ptr = _lib.ebc_be_last_error_msg(self.__be)
        if error_ptr != _ffi.NULL:
            return _ffi.string(error_ptr).decode("utf8")
        else:
            return None

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
            raise PAHMMError("Could not apply GTR model.", self)

        # Tell the garbage collector how to free the resources
        sequences = _ffi.gc(sequences, _lib.ebc_seq_free)

        return Sequences(sequences, self)

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
            raise PAHMMError("Could not apply HKY85 model.", self)

        # Tell the garbage collector how to free the resources
        sequences = _ffi.gc(sequences, _lib.ebc_seq_free)

        return Sequences(sequences, self)

    def execute_jtt_model(self) -> Union[None, "Sequences"]:
        """Generate a Sequences-object that will use the JTT model for distance calculations.

        This model can only be used for amino-acids.
        """

        sequences = _lib.ebc_be_execute_jtt_model(self.__be)

        if sequences == _ffi.NULL:
            raise PAHMMError("Could not apply JTT model.", self)

        # Tell the garbage collector how to free the resources
        sequences = _ffi.gc(sequences, _lib.ebc_seq_free)

        return Sequences(sequences, self)

    def execute_lg_model(self) -> Union[None, "Sequences"]:
        """Generate a Sequences-object that will use the LG model for distance calculations.

        This model can only be used for amino-acids.
        """

        sequences = _lib.ebc_be_execute_lg_model(self.__be)

        if sequences == _ffi.NULL:
            raise PAHMMError("Could not apply LG model.", self)

        # Tell the garbage collector how to free the resources
        sequences = _ffi.gc(sequences, _lib.ebc_seq_free)

        return Sequences(sequences, self)

    def execute_wag_model(self) -> Union[None, "Sequences"]:
        """Generate a Sequences-object that will use the WAG model for distance calculations.

        This model can only be used for amino-acids.
        """

        sequences = _lib.ebc_be_execute_wag_model(self.__be)

        if sequences == _ffi.NULL:
            raise PAHMMError("Could not apply WAG model.", self)

        # Tell the garbage collector how to free the resources
        sequences = _ffi.gc(sequences, _lib.ebc_seq_free)

        return Sequences(sequences, self)

    def set_indel_parameters(self, nb_probability: Union[None, float] = None,
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
        """Decide to load input from a string/bytes/bytearray.

        For better performance use set_file_input instead.
        """

        result: bool = \
            _lib.ebc_be_set_input(self.__be, fasta.encode("utf8") if isinstance(fasta, str) else fasta)

        if not result:
            raise PAHMMError("Could not read FASTA input.", self)

    def set_file_input(self, fasta_path: Union[Path, AnyStr]):
        """Decide to load input from a file.

        This is faster than set_str_input because all data will be read by libpahmm directly,
        bypassing Python.
        """

        fasta_path = self._path_to_bytes(fasta_path)
        result: bool = _lib.ebc_be_set_input_from_file(self.__be, fasta_path)

        if not result:
            raise PAHMMError(f"Could not read sequences from {fasta_path}.", self)


class Sequences:
    """A class used to access sequences and calculate distances between them.

    Note: Do not instantiate an object of this class directly. Use one of
    BandingEstimator's execute_*-methods for that.
    """

    def __init__(self, c_seq, be):
        self._be = be
        self._seq_count = _lib.ebc_seq_count(c_seq)
        self.__seq = c_seq

    def __len__(self):
        """Get the number of sequences.
        """
        return self._seq_count

    def __getitem__(self, seq_id):
        """Get a sequence using a sequence number/ID.
        """
        c_sequence = _lib.ebc_seq_get_sequence(self.__seq, seq_id)
        if c_sequence != _ffi.NULL:
            sequence = _ffi.string(c_sequence)
        else:
            sequence = None

        if self._be.has_last_error():
            raise PAHMMError("Could not get sequence.", self._be)

        return sequence

    def get_distance(self, seq_id1: int, seq_id2: int):
        """Retrieve a distance between two sequences using their numbers/IDs.
        """

        distance = _lib.ebc_seq_get_distance(self.__seq, seq_id1, seq_id2)

        if self._be.has_last_error():
            raise PAHMMError("Could not get distance.", self._be)

        return distance

    def get_distance_from_names(self, seq_name1: Union[bytes, bytearray], seq_name2: Union[bytes, bytearray]):
        """Retrieve a distance between two sequences using their names.
        """
        distance = _lib.ebc_seq_get_distance_from_names(self.__seq, seq_name1, seq_name2)

        if self._be.has_last_error():
            raise PAHMMError("Could not get distance from sequence names.", self._be)

        return distance

    def get_seq_name(self, seq_id: int):
        """Get a sequence using its number or ID.
        """
        c_name = _lib.ebc_seq_get_name(self.__seq, seq_id)
        if c_name != _ffi.NULL:
            name = _ffi.string(c_name)
        else:
            name = None

        if self._be.has_last_error():
            raise PAHMMError("Could not get sequence name.", self._be)

        return name

    def get_sequence_from_name(self, seq_name: Union[bytes, bytearray]):
        """Get a sequence using its name.
        """
        c_name = _lib.ebc_seq_get_sequence_from_name(self.__seq, seq_name)
        if c_name != _ffi.NULL:
            name = _ffi.string(c_name)
        else:
            name = None

        if self._be.has_last_error():
            raise PAHMMError("Could not get sequence from name.", self._be)

        return name
