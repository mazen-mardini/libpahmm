//==============================================================================
// libpahmm - library for paHMM-Tree, a phylogenetic tree estimator
//
// Copyright (c) 2020 Mazen Mardini.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
//     This program is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//     GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses>.
//==============================================================================

#ifndef PAHMM_H
#define PAHMM_H

#include <stdint.h>

#if defined(_MSC_VER) || defined(WIN64) || defined(_WIN64) || defined(__WIN64__) || defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
#  define DECL_EXPORT __declspec(dllexport)
#  define DECL_IMPORT __declspec(dllimport)
#else
#  define DECL_EXPORT     __attribute__((visibility("default")))
#  define DECL_IMPORT     __attribute__((visibility("default")))
#endif

#if defined(PAHMM_LIBRARY)
#  define PAHMM_EXPORT DECL_EXPORT
#else
#  define PAHMM_EXPORT DECL_IMPORT
#endif

#ifndef __cplusplus
#   include <stdbool.h>
#endif

#ifdef __cplusplus
extern "C"
{
#endif

#define EBC_BE_DEFAULTS_INDEL_NB_PROBABILITY 0.0
#define EBC_BE_DEFAULTS_INDEL_RATE 0.0
#define EBC_BE_DEFAULTS_ALPHA 0.5
#define EBC_BE_DEFAULTS_GAMMA_RATE_CATEGORIES 4

    /*
     * The banding estimator used to load sequences from a string
     * or a file and create EBCSequences-objects.
     */
    typedef struct PAHMM_EXPORT EBCBandingEstimator {
        // Parsed input
        void *_parser;
        // Error
        void *_error;

        // Indel params:
        double indel_NB_probability;
        double indel_rate;

        // Discrete Gamma shape parameter alpha:
        double alpha;

        // Gamma rate categories
        unsigned int gamma_rate_categories;

        // Estimate values instead of using predefined parameters
        bool estimate_indel_params;
        bool estimate_alpha;
        bool estimate_categories;
    } EBCBandingEstimator;

    /*
     * A set of sequences.
     */
    typedef struct PAHMM_EXPORT EBCSequences {
        void *_sequences;
        void *_modelEstimator;
        void *_bandingEstimator;
        EBCBandingEstimator *_ebcBandingEstimator;
        int sequenceType;
    } EBCSequences;

    /*
     * Construct a banding estimator object.
     *
     * It is your responsibility to clean it up using ebc_be_free().
     */
    PAHMM_EXPORT EBCBandingEstimator *ebc_be_create();

    /*
     * Destroy a banding estimator object.
     */
    PAHMM_EXPORT void ebc_be_free(EBCBandingEstimator *be);

    /*
     * Get the last error message.
     *
     * If the previous execution was successful, it will return nullptr.
     * The return value will be the error message. You should not free it,
     * it's freed automatically.
     */
    PAHMM_EXPORT const char *ebc_be_last_error_msg(EBCBandingEstimator *be);

    /*
     * Use the ebc_be_execute_* functions to actually perform distance
     * calculations between the input sequences.
     *
     * The return value is an EBCSequences-object upon successful execution,
     * and it is your responsibility to clean them up using ebc_seq_free().
     * If an error occurred, it will return NULL.
     */

    // General time reversible substitution (GTR) model
    PAHMM_EXPORT EBCSequences *ebc_be_execute_gtr_model(EBCBandingEstimator *be,
                                    double, double, double, double, double);
    PAHMM_EXPORT EBCSequences *ebc_be_execute_gtr_modelv2(EBCBandingEstimator *be);

    // HKY85 substitution model
    PAHMM_EXPORT EBCSequences *ebc_be_execute_hky85_model(EBCBandingEstimator *be, double);
    PAHMM_EXPORT EBCSequences *ebc_be_execute_hky85_modelv2(EBCBandingEstimator *be);

    // Jones 1992 AA substitution model
    PAHMM_EXPORT EBCSequences *ebc_be_execute_jtt_model(EBCBandingEstimator *be);

    // Le & Gasquel AA substitution model
    PAHMM_EXPORT EBCSequences *ebc_be_execute_lg_model(EBCBandingEstimator *be);

    // Whelan and Goldman AA substitution model
    PAHMM_EXPORT EBCSequences *ebc_be_execute_wag_model(EBCBandingEstimator *be);

    /*
     * Set/Unset different parameters
     *
     * Use ebc_be_set_* or ebc_be_unset_* functions
     * to set parameters, or unset parameters and let the model estimator
     * choose their values.
     */

    // Indel parameters
    PAHMM_EXPORT void ebc_be_set_indel_parameters(EBCBandingEstimator *be,
                                     double NB_probability, double rate);
    PAHMM_EXPORT void ebc_be_unset_indel_parameters(EBCBandingEstimator *be);

    // Discrete Gamma shape parameter alpha
    PAHMM_EXPORT void ebc_be_set_alpha(EBCBandingEstimator *be, double alpha);
    PAHMM_EXPORT void ebc_be_unset_alpha(EBCBandingEstimator *be);

    // Gamma rate categories
    PAHMM_EXPORT void ebc_be_set_categories(EBCBandingEstimator *be, unsigned int categories);
    PAHMM_EXPORT void ebc_be_unset_categories(EBCBandingEstimator *be);

    /*
     * Set sequence input. Should be in FASTA-format.
     *
     * These functions return false upon failure, and true when they are successful.
     */

    // Set sequence input using a string
    PAHMM_EXPORT bool ebc_be_set_input(EBCBandingEstimator *be, const char *fasta);

    // Set sequence input from a file
    PAHMM_EXPORT bool ebc_be_set_input_from_file(EBCBandingEstimator *be, const char *file_name);

    /*
     * Destroy a sequences object.
     */
    PAHMM_EXPORT void ebc_seq_free(EBCSequences *seq);

    /*
     * Get sequence count.
     *
     * If, for example count = 5, then all existing sequence ID's are 0, 1, 2, 3 and 4.
     */
    PAHMM_EXPORT unsigned int ebc_seq_count(EBCSequences *seq);

    /*
     * Get the distance between two sequences.
     *
     * If an error occurs, NAN is returned.
     *
     * If the distance hasn't been calulated before, it will be calculated and the function
     * will return the result.
     */
    PAHMM_EXPORT double ebc_seq_get_distance(EBCSequences *seq, unsigned int seq_id1, unsigned int seq_id2);

    /*
     * Get the distance between two sequences using their names.
     *
     * If an error occurs, NAN is returned.
     *
     * If the distance hasn't been calulated before, it will be calculated and the function
     * will return the result.
     */
    PAHMM_EXPORT double ebc_seq_get_distance_from_names(EBCSequences *seq,
                                                        const char *seq_name1, const char *seq_name2);

    /*
     * Get the name of a sequence from a sequence ID.
     *
     * If an error occurs, NULL is returned.
     *
     * Note: You do not have to worry about freeing the return value, this is done by
     * ebc_seq_free().
     */
    PAHMM_EXPORT const char *ebc_seq_get_name(EBCSequences *seq, unsigned int seq_id);

    /*
     * Get the sequence string from a sequence ID.
     *
     * If an error occurs, NULL is returned.
     *
     * Note: You do not have to worry about freeing the return value, this is done by
     * ebc_seq_free().
     */
    PAHMM_EXPORT const char *ebc_seq_get_sequence(EBCSequences *seq, unsigned int seq_id);

    /*
     * Get the sequence string from a sequence name.
     *
     * If an error occurs, NULL is returned.
     *
     * Note: You do not have to worry about freeing the return value, this is done by
     * ebc_seq_free().
     */
    PAHMM_EXPORT const char *ebc_seq_get_sequence_from_name(EBCSequences *seq, const char *seq_name);


#ifdef __cplusplus
}
#endif

#endif // PAHMM_H
