/*************************************************************************************
                           The MIT License

   BWA-MEM2  (Sequence alignment using Burrows-Wheeler Transform),
   Copyright (C) 2019  Intel Corporation, Heng Li.

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.

Authors: Sanchit Misra <sanchit.misra@intel.com>; Vasimuddin Md <vasimuddin.md@intel.com>;
*****************************************************************************************/

// submodule test comment

#ifndef _FMI_SEARCH_H
#define _FMI_SEARCH_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <immintrin.h>
#include <limits.h>
#include <fstream>

#include "read_index_ele.h"
#include "bwa.h"

#define DUMMY_CHAR 6

#define assert_not_null(x, size, cur_alloc) \
        if (x == NULL) { fprintf(stderr, "Allocation of %0.2lf GB for " #x " failed.\nCurrent Allocation = %0.2lf GB\n", size * 1.0 /(1024*1024*1024), cur_alloc * 1.0 /(1024*1024*1024)); exit(EXIT_FAILURE); }

// TODO: create new variable OCC_COMPRESSION_FACTOR and generate relavant vars accordingly
#define CP_BLOCK_SIZE 1
#define CP_FILENAME_SUFFIX ".bwt.2bit.1"
#define CP_MASK 0
#define CP_SHIFT 0

//#define CP_BLOCK_SIZE 8
//#define CP_FILENAME_SUFFIX ".bwt.2bit.8"
//#define CP_MASK 7
//#define CP_SHIFT 3

//#define CP_BLOCK_SIZE 16
//#define CP_FILENAME_SUFFIX ".bwt.2bit.16"
//#define CP_MASK 15
//#define CP_SHIFT 4

//#define CP_BLOCK_SIZE 32
//#define CP_FILENAME_SUFFIX ".bwt.2bit.32"
//#define CP_MASK 31
//#define CP_SHIFT 5

//#define CP_BLOCK_SIZE 64
//#define CP_FILENAME_SUFFIX ".bwt.2bit.64"
//#define CP_MASK 63
//#define CP_SHIFT 6

typedef struct checkpoint_occ_scalar_uncompressed
{
    int64_t cp_count[4];
    uint8_t bwt_char;
}CP_OCC_UNCOMPRESSED;

// GET_OCC:
// pp is the occ table row that we wish to retrieve, c is the base.
// occ_id_pp, y_pp, one_hot_bwt_str_c_pp, and match_mask_pp, are temporary variables used only in GET_OCC, 
// occ_pp is the "return" value of this macro, i.e. the occurance value retrieved from the occ table.
//
// Example: pp = 1 = 0b001, compression ration = 4, CP_SHIFT = 2, CP_MASK = 3 = 0b11
// This means that every 4th entry is stored in the occ table: 0th, 4th, 8th, 12th rows etc.

// e.g. uncompressed Occ table
// ----------------------
// row | bwt | count[C] |
// 0   |  A  | 22       |
// 1   |  C  | 23       |
// 2   |  C  | 24       |
// 3   |  G  | 24       |
// 4   |  A  | 24       |
// 5   |  G  | 24       |
// 6   |  T  | 24       |
// 7   |  T  | 24       |
// ----------------------
// The corresponding compressed Occ table, compression factor = 4
// --------------------------------------------------------------------
// row | count[C] | bwt substring | one_hot_bwt_str[A] | C  | G  | T  |
// 0   |   22     |  ACCG         |       1000         |0110|0001|0000|
// 1   |   24     |  AGTT         |       1000         |0000|0100|0011|             
// --------------------------------------------------------------------

// occ_id_pp = 0b0 = 0, ypp = 0b10 = 1
// occ_pp = cp_occ[1].cp_count = the 0th row of uncompressed table 
// Since we are looking for row 1, need to used count[C] at 0th row to compute count[C] at row 1
// one_hot_bwt_str_c_pp = cp_occ[0].one_hot_bwt_str[C] = 0110
// See FM_index.cpp one_hot_bwt_str construction to see how one_hot_bwt_str works. Also see BWA-MEM2 paper section IV.A
// match_mask_pp = 0110 & one_hot_mask_array[1]
// This mask is used to only count the number of Cs up to row 1. Thus, one_hot_mask_array should be 1100.
// match_mask_pp = 0100
// _mm_countbits_ counts the number of 1 bits in match_mask_pp
// occ_pp = 22 + 1 = 23 --> occ[C, 1] = 23


// If no compression:
// Just return cp_count of the requested row
// No need to store one_hot_bwt_str[] for each row

// TODO: dont need occ_id_pp anymore
#define \
GET_OCC_UNCOMPRESSED(pp, c, occ_id_pp, occ_pp) \
                int64_t occ_id_pp = pp; \
                int64_t occ_pp = cp_occ[pp].cp_count[c];
                

typedef struct smem_struct
{
#ifdef DEBUG
    uint64_t info; // for debug
#endif
    uint32_t rid;
    uint32_t m, n;
    int64_t k, l, s;
}SMEM;

#define SAL_PFD 16

class FMI_search: public indexEle
{
    public:
    FMI_search(const char *fname);
    ~FMI_search();
    //int64_t beCalls;
    
    int build_index();
    void load_index();

    void getSMEMs(uint8_t *enc_qdb,
                  int32_t numReads,
                  int32_t batch_size,
                  int32_t readlength,
                  int32_t minSeedLengh,
                  int32_t numthreads,
                  SMEM *matchArray,
                  int64_t *numTotalSmem);
    
    void getSMEMsOnePosOneThread(uint8_t *enc_qdb,
                                 int16_t *query_pos_array,
                                 int32_t *min_intv_array,
                                 int32_t *rid_array,
                                 int32_t numReads,
                                 int32_t batch_size,
                                 const bseq1_t *seq_,
                                 int32_t *query_cum_len_ar,
                                 int32_t  max_readlength,
                                 int32_t minSeedLen,
                                 SMEM *matchArray,
                                 int64_t *__numTotalSmem);
    
    void getSMEMsAllPosOneThread(uint8_t *enc_qdb,
                                 int32_t *min_intv_array,
                                 int32_t *rid_array,
                                 int32_t numReads,
                                 int32_t batch_size,
                                 const bseq1_t *seq_,
                                 int32_t *query_cum_len_ar,
                                 int32_t max_readlength,
                                 int32_t minSeedLen,
                                 SMEM *matchArray,
                                 int64_t *__numTotalSmem);
        
    
    int64_t bwtSeedStrategyAllPosOneThread(uint8_t *enc_qdb,
                                           int32_t *max_intv_array,
                                           int32_t numReads,
                                           const bseq1_t *seq_,
                                           int32_t *query_cum_len_ar,
                                           int32_t minSeedLen,
                                           SMEM *matchArray);
        
    void sortSMEMs(SMEM *matchArray,
                   int64_t numTotalSmem[],
                   int32_t numReads,
                   int32_t readlength,
                   int nthreads);
    int64_t get_sa_entry(int64_t pos);
    void get_sa_entries(int64_t *posArray,
                        int64_t *coordArray,
                        uint32_t count,
                        int32_t nthreads);
    void get_sa_entries(SMEM *smemArray,
                        int64_t *coordArray,
                        int32_t *coordCountArray,
                        uint32_t count,
                        int32_t max_occ);
    int64_t get_sa_entry_compressed(int64_t pos, int tid=0);
    void get_sa_entries(SMEM *smemArray,
                        int64_t *coordArray,
                        int32_t *coordCountArray,
                        uint32_t count,
                        int32_t max_occ,
                        int tid);
    int64_t call_one_step(int64_t pos, int64_t &sa_entry, int64_t &offset);
    void get_sa_entries_prefetch(SMEM *smemArray, int64_t *coordArray,
                                 int64_t *coordCountArray, int64_t count,
                                 const int32_t max_occ, int tid, int64_t &id_);
    
    int64_t reference_seq_len;
    int64_t sentinel_index;
private:
        char file_name[PATH_MAX];
        int64_t index_alloc;
        int64_t count[5];
        uint32_t *sa_ls_word;
        int8_t *sa_ms_byte;
        CP_OCC_UNCOMPRESSED *cp_occ;

        int64_t pac_seq_len(const char *fn_pac);
        void pac2nt(const char *fn_pac,
                    std::string &reference_seq);
        int build_fm_index(const char *ref_file_name,
                               char *binary_seq,
                               int64_t ref_seq_len,
                               int64_t *sa_bwt,
                               int64_t *count);
        SMEM backwardExt(SMEM smem, uint8_t a);
};

#endif
