from cffi import FFI
from python.paths import *
import sys

ffibuilder = FFI()
ffibuilder.cdef("""
    typedef struct EBCSequences {
        void *_sequences;
        void *_modelEstimator;
        void *_bandingEstimator;
        int sequenceType;
    } EBCSequences;
    typedef struct EBCBandingEstimator {
        void *_parser;
        double indel_NB_probability;
        double indel_rate;
        double alpha;
        unsigned int gamma_rate_categories;
        bool estimate_indel_params;
        bool estimate_alpha;
        bool estimate_categories;
    } EBCBandingEstimator;
    EBCBandingEstimator *ebc_be_create();
    void ebc_be_free(EBCBandingEstimator *);
    const char *ebc_be_last_error_msg(EBCBandingEstimator *be);
    EBCSequences *ebc_be_execute_gtr_model(EBCBandingEstimator *be,
                                    double, double, double, double, double);
    EBCSequences *ebc_be_execute_gtr_modelv2(EBCBandingEstimator *be);
    EBCSequences *ebc_be_execute_hky85_model(EBCBandingEstimator *be, double);
    EBCSequences *ebc_be_execute_hky85_modelv2(EBCBandingEstimator *be);
    EBCSequences *ebc_be_execute_jtt_model(EBCBandingEstimator *be);
    EBCSequences *ebc_be_execute_lg_model(EBCBandingEstimator *be);
    EBCSequences *ebc_be_execute_wag_model(EBCBandingEstimator *be);
    void ebc_be_set_indel_parameters(EBCBandingEstimator *be,
                                     double NB_probability, double rate);
    void ebc_be_unset_indel_parameters(EBCBandingEstimator *be);
    void ebc_be_set_alpha(EBCBandingEstimator *be, double alpha);
    void ebc_be_unset_alpha(EBCBandingEstimator *be);
    void ebc_be_set_categories(EBCBandingEstimator *be, unsigned int categories);
    void ebc_be_unset_categories(EBCBandingEstimator *be);
    bool ebc_be_set_input(EBCBandingEstimator *be, const char *fasta);
    bool ebc_be_set_input_from_file(EBCBandingEstimator *be, const char *file_name);
    void ebc_seq_free(EBCSequences *seq);
    unsigned int ebc_seq_count(EBCSequences *seq);
    double ebc_seq_get_distance(EBCSequences *seq, unsigned int seq_id1, unsigned int seq_id2);
    double ebc_seq_get_distance_from_names(EBCSequences *seq, const char *seq_name1, const char *seq_name2);
    const char *ebc_seq_get_name(EBCSequences *seq, unsigned int seq_id);
    const char *ebc_seq_get_sequence(EBCSequences *seq, unsigned int seq_id);
    const char *ebc_seq_get_sequence_from_name(EBCSequences *seq, const char *seq_name);
""")

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
