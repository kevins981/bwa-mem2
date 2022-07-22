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

#include <stdio.h>
#include "sais.h"
#include "FMI_search.h"
#include "memcpy_bwamem.h"
#include "profiling.h"

#ifdef __cplusplus
extern "C" {
#endif
#include "safe_str_lib.h"
#ifdef __cplusplus
}
#endif

FMI_search::FMI_search(const char *fname)
{
    fprintf(stderr, "* Entering FMI_search\n");
    //strcpy(file_name, fname);
    strcpy_s(file_name, PATH_MAX, fname);
    reference_seq_len = 0;
    sentinel_index = 0;
    index_alloc = 0;
    sa_ls_word = NULL;
    sa_ms_byte = NULL;
    cp_occ = NULL;
}

FMI_search::~FMI_search()
{
    if(sa_ms_byte)
        _mm_free(sa_ms_byte);
    if(sa_ls_word)
        _mm_free(sa_ls_word);
    if(cp_occ)
        _mm_free(cp_occ);
}

int64_t FMI_search::pac_seq_len(const char *fn_pac)
{
	FILE *fp;
	int64_t pac_len;
	uint8_t c;
	fp = xopen(fn_pac, "rb");
	err_fseek(fp, -1, SEEK_END);
	pac_len = err_ftell(fp);
	err_fread_noeof(&c, 1, 1, fp);
	err_fclose(fp);
	return (pac_len - 1) * 4 + (int)c;
}

void FMI_search::pac2nt(const char *fn_pac, std::string &reference_seq)
{
	uint8_t *buf2;
	int64_t i, pac_size, seq_len;
	FILE *fp;

	// initialization
	seq_len = pac_seq_len(fn_pac);
    assert(seq_len > 0);
    assert(seq_len <= 0x7fffffffffL);
	fp = xopen(fn_pac, "rb");

	// prepare sequence
	pac_size = (seq_len>>2) + ((seq_len&3) == 0? 0 : 1);
	buf2 = (uint8_t*)calloc(pac_size, 1);
    assert(buf2 != NULL);
	err_fread_noeof(buf2, 1, pac_size, fp);
	err_fclose(fp);
	for (i = 0; i < seq_len; ++i) {
		int nt = buf2[i>>2] >> ((3 - (i&3)) << 1) & 3;
        switch(nt)
        {
            case 0:
                reference_seq += "A";
            break;
            case 1:
                reference_seq += "C";
            break;
            case 2:
                reference_seq += "G";
            break;
            case 3:
                reference_seq += "T";
            break;
            default:
                fprintf(stderr, "ERROR! Value of nt is not in 0,1,2,3!");
                exit(EXIT_FAILURE);
        }
	}
    for(i = seq_len - 1; i >= 0; i--)
    {
        char c = reference_seq[i];
        switch(c)
        {
            case 'A':
                reference_seq += "T";
            break;
            case 'C':
                reference_seq += "G";
            break;
            case 'G':
                reference_seq += "C";
            break;
            case 'T':
                reference_seq += "A";
            break;
        }
    }
	free(buf2);
}

int FMI_search::build_fm_index(const char *ref_file_name, char *binary_seq, int64_t ref_seq_len, int64_t *sa_bwt, int64_t *count) {
    printf("ref_seq_len = %ld\n", ref_seq_len);
    fflush(stdout);

    char outname[PATH_MAX];

    strcpy_s(outname, PATH_MAX, ref_file_name);
    strcat_s(outname, PATH_MAX, CP_FILENAME_SUFFIX);
    //sprintf(outname, "%s.bwt.2bit.%d", ref_file_name, CP_BLOCK_SIZE);

    std::fstream outstream (outname, std::ios::out | std::ios::binary);
    outstream.seekg(0);	

    printf("count = %ld, %ld, %ld, %ld, %ld\n", count[0], count[1], count[2], count[3], count[4]);
    fflush(stdout);

    uint8_t *bwt;

    ref_seq_len++;
    outstream.write((char *)(&ref_seq_len), 1 * sizeof(int64_t));
    outstream.write((char*)count, 5 * sizeof(int64_t));

    int64_t i;
    int64_t ref_seq_len_aligned = ((ref_seq_len + CP_BLOCK_SIZE - 1) / CP_BLOCK_SIZE) * CP_BLOCK_SIZE;
    int64_t size = ref_seq_len_aligned * sizeof(uint8_t);
    bwt = (uint8_t *)_mm_malloc(size, 64);
    assert_not_null(bwt, size, index_alloc);

    int64_t sentinel_index = -1;
    for(i=0; i< ref_seq_len; i++)
    {
        if(sa_bwt[i] == 0)
        {
            bwt[i] = 4;
            printf("BWT[%ld] = 4\n", i);
            sentinel_index = i;
        }
        else
        {
            char c = binary_seq[sa_bwt[i]-1];
            switch(c)
            {
                case 0: bwt[i] = 0;
                          break;
                case 1: bwt[i] = 1;
                          break;
                case 2: bwt[i] = 2;
                          break;
                case 3: bwt[i] = 3;
                          break;
                default:
                        fprintf(stderr, "ERROR! i = %ld, c = %c\n", i, c);
                        exit(EXIT_FAILURE);
            }
        }
    }
    for(i = ref_seq_len; i < ref_seq_len_aligned; i++)
        bwt[i] = DUMMY_CHAR;


    printf("CP_SHIFT = %d, CP_MASK = %d\n", CP_SHIFT, CP_MASK);
    printf("sizeof CP_OCC_UNCOMPRESSED = %ld\n", sizeof(CP_OCC_UNCOMPRESSED));
    fflush(stdout);
    // create checkpointed occ
    int64_t cp_occ_size = (ref_seq_len >> CP_SHIFT) + 1;
    CP_OCC_UNCOMPRESSED *cp_occ = NULL;

    size = cp_occ_size * sizeof(CP_OCC_UNCOMPRESSED);
    cp_occ = (CP_OCC_UNCOMPRESSED *)_mm_malloc(size, 64);
    assert_not_null(cp_occ, size, index_alloc);
    memset(cp_occ, 0, cp_occ_size * sizeof(CP_OCC_UNCOMPRESSED));
    int64_t cp_count[16];

    memset(cp_count, 0, 16 * sizeof(int64_t));
    // TODO: upper bound should be cp_occ_size, which is ref_seq_len+1
    for(i = 0; i < ref_seq_len; i++)
    {
        // Only store Occ table rows that are multiples of the compression factor
        CP_OCC_UNCOMPRESSED cpo; // one Occ table row
        cpo.cp_count[0] = cp_count[0];
        cpo.cp_count[1] = cp_count[1];
        cpo.cp_count[2] = cp_count[2];
        cpo.cp_count[3] = cp_count[3];

        cpo.bwt_char = bwt[i];

        cp_occ[i >> CP_SHIFT] = cpo;

        cp_count[bwt[i]]++;
    }
    outstream.write((char*)cp_occ, cp_occ_size * sizeof(CP_OCC_UNCOMPRESSED));
    _mm_free(cp_occ);
    _mm_free(bwt);

    #if SA_COMPRESSION  

    size = ((ref_seq_len >> SA_COMPX)+ 1)  * sizeof(uint32_t);
    uint32_t *sa_ls_word = (uint32_t *)_mm_malloc(size, 64);
    assert_not_null(sa_ls_word, size, index_alloc);
    size = ((ref_seq_len >> SA_COMPX) + 1) * sizeof(int8_t);
    int8_t *sa_ms_byte = (int8_t *)_mm_malloc(size, 64);
    assert_not_null(sa_ms_byte, size, index_alloc);
    int64_t pos = 0;
    // Building compressed suffix array. The SA is the same length as the reference sequence
    for(i = 0; i < ref_seq_len; i++)
    {
        if ((i & SA_COMPX_MASK) == 0)
        {
            // If SA index is a multiple of compression factor, store in memory. Else, do not store.
            // This is the same compression scheme as the occ table. 
            sa_ls_word[pos] = sa_bwt[i] & 0xffffffff; // lower 32 bits
            sa_ms_byte[pos] = (sa_bwt[i] >> 32) & 0xff; // upper 8 bits
            pos++;
        }
    }
    fprintf(stderr, "pos: %d, ref_seq_len__: %ld\n", pos, ref_seq_len >> SA_COMPX);
    outstream.write((char*)sa_ms_byte, ((ref_seq_len >> SA_COMPX) + 1) * sizeof(int8_t));
    outstream.write((char*)sa_ls_word, ((ref_seq_len >> SA_COMPX) + 1) * sizeof(uint32_t));
    
    #else
    
    size = ref_seq_len * sizeof(uint32_t);
    uint32_t *sa_ls_word = (uint32_t *)_mm_malloc(size, 64);
    assert_not_null(sa_ls_word, size, index_alloc);
    size = ref_seq_len * sizeof(int8_t);
    int8_t *sa_ms_byte = (int8_t *)_mm_malloc(size, 64);
    assert_not_null(sa_ms_byte, size, index_alloc);
    for(i = 0; i < ref_seq_len; i++)
    {
        sa_ls_word[i] = sa_bwt[i] & 0xffffffff;
        sa_ms_byte[i] = (sa_bwt[i] >> 32) & 0xff;
    }
    outstream.write((char*)sa_ms_byte, ref_seq_len * sizeof(int8_t));
    outstream.write((char*)sa_ls_word, ref_seq_len * sizeof(uint32_t));
    
    #endif

    outstream.write((char *)(&sentinel_index), 1 * sizeof(int64_t));
    outstream.close();
    printf("max_occ_ind = %ld\n", i >> CP_SHIFT);    
    fflush(stdout);

    _mm_free(sa_ms_byte);
    _mm_free(sa_ls_word);
    return 0;
}

int FMI_search::build_index() {
    printf("CP_SHIFT = %d, CP_MASK = %d\n", CP_SHIFT, CP_MASK);

    char *prefix = file_name;
    uint64_t startTick;
    startTick = __rdtsc();
    index_alloc = 0;

    std::string reference_seq;
    char pac_file_name[PATH_MAX];
    strcpy_s(pac_file_name, PATH_MAX, prefix);
    strcat_s(pac_file_name, PATH_MAX, ".pac");
    //sprintf(pac_file_name, "%s.pac", prefix);
    pac2nt(pac_file_name, reference_seq);
	int64_t pac_len = reference_seq.length();
    int status;
    int64_t size = pac_len * sizeof(char);
    char *binary_ref_seq = (char *)_mm_malloc(size, 64);
    index_alloc += size;
    assert_not_null(binary_ref_seq, size, index_alloc);
    char binary_ref_name[PATH_MAX];
    strcpy_s(binary_ref_name, PATH_MAX, prefix);
    strcat_s(binary_ref_name, PATH_MAX, ".0123");
    //sprintf(binary_ref_name, "%s.0123", prefix);
    std::fstream binary_ref_stream (binary_ref_name, std::ios::out | std::ios::binary);
    binary_ref_stream.seekg(0);
    fprintf(stderr, "init ticks = %llu\n", __rdtsc() - startTick);
    startTick = __rdtsc();
    int64_t i, count[16];
	memset(count, 0, sizeof(int64_t) * 16);
    for(i = 0; i < pac_len; i++)
    {
        switch(reference_seq[i])
        {
            case 'A':
            binary_ref_seq[i] = 0, ++count[0];
            break;
            case 'C':
            binary_ref_seq[i] = 1, ++count[1];
            break;
            case 'G':
            binary_ref_seq[i] = 2, ++count[2];
            break;
            case 'T':
            binary_ref_seq[i] = 3, ++count[3];
            break;
            default:
            binary_ref_seq[i] = 4;

        }
    }
    count[4]=count[0]+count[1]+count[2]+count[3];
    count[3]=count[0]+count[1]+count[2];
    count[2]=count[0]+count[1];
    count[1]=count[0];
    count[0]=0;
    fprintf(stderr, "ref seq len = %ld\n", pac_len);
    binary_ref_stream.write(binary_ref_seq, pac_len * sizeof(char));
    fprintf(stderr, "binary seq ticks = %llu\n", __rdtsc() - startTick);
    startTick = __rdtsc();

    size = (pac_len + 2) * sizeof(int64_t);
    int64_t *suffix_array=(int64_t *)_mm_malloc(size, 64);
    index_alloc += size;
    assert_not_null(suffix_array, size, index_alloc);
    startTick = __rdtsc();
	//status = saisxx<const char *, int64_t *, int64_t>(reference_seq.c_str(), suffix_array + 1, pac_len, 4);
	status = saisxx(reference_seq.c_str(), suffix_array + 1, pac_len);
	suffix_array[0] = pac_len;
    fprintf(stderr, "build suffix-array ticks = %llu\n", __rdtsc() - startTick);
    startTick = __rdtsc();

	build_fm_index(prefix, binary_ref_seq, pac_len, suffix_array, count);
    fprintf(stderr, "build fm-index ticks = %llu\n", __rdtsc() - startTick);
    _mm_free(binary_ref_seq);
    _mm_free(suffix_array);
    return 0;
}

void FMI_search::load_index()
{
    char *ref_file_name = file_name;
    //beCalls = 0;
    char cp_file_name[PATH_MAX];
    strcpy_s(cp_file_name, PATH_MAX, ref_file_name);
    // cp_file_name is broad.bwt.2bit.64A. This sounds like the BWT vector,
    // i.e. last column of the BW matrix
    strcat_s(cp_file_name, PATH_MAX, CP_FILENAME_SUFFIX);

    // Read the BWT and FM index of the reference sequence
    // I am not sure what the above comment means. The BWT is a part of the FM index.
    FILE *cpstream = NULL;
    cpstream = fopen(cp_file_name,"rb");
    fprintf(stderr, "Opening cp file: %s\n", cp_file_name);

    if (cpstream == NULL)
    {
        fprintf(stderr, "ERROR! Unable to open the file: %s\n", cp_file_name);
        exit(EXIT_FAILURE);
    }
    else
    {
        fprintf(stderr, "* Index file found. Loading index from %s\n", cp_file_name);
    }

    // fread(ptr, size, cnt, stream); read size*cnt into ptr
    err_fread_noeof(&reference_seq_len, sizeof(int64_t), 1, cpstream);
    // This file contains the reference sequence length. 
    assert(reference_seq_len > 0);
    assert(reference_seq_len <= 0x7fffffffffL);

    // For large, this is 6434693835, or 6.4G
    fprintf(stderr, "* Reference seq len for bi-index = %ld\n", reference_seq_len);

    // create checkpointed occ
    // "checkpoint" in this context means precomputing the tally.
    int64_t cp_occ_size = (reference_seq_len >> CP_SHIFT) + 1;
    // CP_SHIFT is 6, so shift right means divide by 2^6 = 64.
    // The full OCC table should have the same size as the reference sequence.
    // The divide by 64 may be because we are only keeping one row for every
    // 64 rows of the full OCC table (in order to reduce memory usage).
    cp_occ = NULL;

    // count is the data structure that describes the first column of the BW matrix
    err_fread_noeof(&count[0], sizeof(int64_t), 5, cpstream);
    if ((cp_occ = (CP_OCC_UNCOMPRESSED *)_mm_malloc(cp_occ_size * sizeof(CP_OCC_UNCOMPRESSED), 64)) == NULL) {
        // cp_occ_size = 100542092 or 100M
        // sizeof(CP_OCC) = 64. Wait how does that work. is this just a pointer?
        // sizeof(CP_OCC_UNCOMPRESSED) = ?. 
        // _mm_malloc(size, 64) allocates 64-byte alighed memory 
        // so malloc 6.4GB of memory
        fprintf(stderr, "ERROR! unable to allocated cp_occ memory\n");
        exit(EXIT_FAILURE);
    }

    // read all of 6.4GB OCC table. Seems like this file contains the OCC
    err_fread_noeof(cp_occ, sizeof(CP_OCC_UNCOMPRESSED), cp_occ_size, cpstream);
    // Workaround for uncompressed: the algorithm is access cp_occ[cp_occ_size-1], which we did not set
    // Adding this last entry manually. Hardcoding for now.
    // TODO: properly generate the last entry in cp_occ during build_index()
    cp_occ[cp_occ_size-1].cp_count[0] = 1882204623;
    cp_occ[cp_occ_size-1].cp_count[1] = 1335142294;
    cp_occ[cp_occ_size-1].cp_count[2] = 1335142294;
    cp_occ[cp_occ_size-1].cp_count[3] = 1882204623;

    int64_t ii = 0;
    for(ii = 0; ii < 5; ii++)// update read count structure
    {
        count[ii] = count[ii] + 1;
    }

    // compression is enabled by default
    #if SA_COMPRESSION

    int64_t reference_seq_len_ = (reference_seq_len >> SA_COMPX) + 1;
    sa_ms_byte = (int8_t *)_mm_malloc(reference_seq_len_ * sizeof(int8_t), 64);
    sa_ls_word = (uint32_t *)_mm_malloc(reference_seq_len_ * sizeof(uint32_t), 64);
    err_fread_noeof(sa_ms_byte, sizeof(int8_t), reference_seq_len_, cpstream);
    err_fread_noeof(sa_ls_word, sizeof(uint32_t), reference_seq_len_, cpstream);
    
    #else
    
    sa_ms_byte = (int8_t *)_mm_malloc(reference_seq_len * sizeof(int8_t), 64);
    sa_ls_word = (uint32_t *)_mm_malloc(reference_seq_len * sizeof(uint32_t), 64);
    err_fread_noeof(sa_ms_byte, sizeof(int8_t), reference_seq_len, cpstream);
    err_fread_noeof(sa_ls_word, sizeof(uint32_t), reference_seq_len, cpstream);

    #endif

    sentinel_index = -1;
    #if SA_COMPRESSION
    err_fread_noeof(&sentinel_index, sizeof(int64_t), 1, cpstream);
    // not sure what this sentinel index is. Sentinel refers to the terminating $ character
    fprintf(stderr, "* sentinel-index: %ld\n", sentinel_index);
    #endif
    fclose(cpstream);

    int64_t x;
    #if !SA_COMPRESSION
    for(x = 0; x < reference_seq_len; x++)
    {
        // fprintf(stderr, "x: %ld\n", x);
        // what is happening here
        #if SA_COMPX_MASK
        if(get_sa_entry_compressed(x) == 0) {
            sentinel_index = x;
            break;
        }
        #else
        if(get_sa_entry(x) == 0) {
            sentinel_index = x;
            break;
        }
        #endif
    }
    fprintf(stderr, "\nsentinel_index: %ld\n", x);    
    #endif

    fprintf(stderr, "* Count:\n");
    for(x = 0; x < 5; x++)
    {
        fprintf(stderr, "%ld,\t%lu\n", x, (unsigned long)count[x]);
    }
    fprintf(stderr, "\n");  
    // count content:
    // 0,	1
    // 1,	1882204624
    // 2,	3217346918
    // 3,	4552489212
    // 4,	6434693835
    // 1234 probably represents ACGT (in a particular order). 
    // So in the BW matrix, row 1 to row 1882204624 is probably one letter.

    fprintf(stderr, "* Reading other elements of the index from files %s\n",
            ref_file_name); // /ssd1/kevin/genomics/input-datasets/fmi/broad
    bwa_idx_load_ele(ref_file_name, BWA_IDX_ALL);

    fprintf(stderr, "* Done reading Index!!\n");
}

void FMI_search::getSMEMsOnePosOneThread(uint8_t *enc_qdb,
                                         int16_t *query_pos_array,
                                         int32_t *min_intv_array,
                                         int32_t *rid_array,
                                         int32_t numReads,
                                         int32_t batch_size,
                                         const bseq1_t *seq_,
                                         int32_t *query_cum_len_ar,
                                         int32_t max_readlength,
                                         int32_t minSeedLen,
                                         SMEM *matchArray,
                                         int64_t *__numTotalSmem)
{
    int64_t numTotalSmem = *__numTotalSmem;
    SMEM prevArray[max_readlength];

    uint32_t i;
    // Perform SMEM for original reads
    for(i = 0; i < numReads; i++)
    {
        int x = query_pos_array[i];
        int32_t rid = rid_array[i];
        int next_x = x + 1;

        int readlength = seq_[rid].l_seq;
        int offset = query_cum_len_ar[rid];
        // uint8_t a = enc_qdb[rid * readlength + x];
        uint8_t a = enc_qdb[offset + x];

        if(a < 4)
        {
            SMEM smem;
            smem.rid = rid;
            smem.m = x;
            smem.n = x;
            smem.k = count[a];
            smem.l = count[3 - a];
            smem.s = count[a+1] - count[a];
            int numPrev = 0;
            
            int j;
            for(j = x + 1; j < readlength; j++)
            {
                // a = enc_qdb[rid * readlength + j];
                a = enc_qdb[offset + j];
                next_x = j + 1;
                if(a < 4)
                {
                    SMEM smem_ = smem;

                    // Forward extension is backward extension with the BWT of reverse complement
                    smem_.k = smem.l;
                    smem_.l = smem.k;
                    SMEM newSmem_ = backwardExt(smem_, 3 - a);
                    //SMEM newSmem_ = forwardExt(smem_, 3 - a);
                    SMEM newSmem = newSmem_;
                    newSmem.k = newSmem_.l;
                    newSmem.l = newSmem_.k;
                    newSmem.n = j;

                    int32_t s_neq_mask = newSmem.s != smem.s;

                    prevArray[numPrev] = smem;
                    numPrev += s_neq_mask;
                    if(newSmem.s < min_intv_array[i])
                    {
                        next_x = j;
                        break;
                    }
                    smem = newSmem;
#ifdef ENABLE_PREFETCH
                    _mm_prefetch((const char *)(&cp_occ[(smem.k) >> CP_SHIFT]), _MM_HINT_T0);
                    _mm_prefetch((const char *)(&cp_occ[(smem.l) >> CP_SHIFT]), _MM_HINT_T0);
#endif
                }
                else
                {
                    break;
                }
            }
            if(smem.s >= min_intv_array[i])
            {

                prevArray[numPrev] = smem;
                numPrev++;
            }

            SMEM *prev;
            prev = prevArray;

            int p;
            for(p = 0; p < (numPrev/2); p++)
            {
                SMEM temp = prev[p];
                prev[p] = prev[numPrev - p - 1];
                prev[numPrev - p - 1] = temp;
            }

            // Backward search
            int cur_j = readlength;
            for(j = x - 1; j >= 0; j--)
            {
                int numCurr = 0;
                int curr_s = -1;
                // a = enc_qdb[rid * readlength + j];
                a = enc_qdb[offset + j];

                if(a > 3)
                {
                    break;
                }
                for(p = 0; p < numPrev; p++)
                {
                    SMEM smem = prev[p];
                    SMEM newSmem = backwardExt(smem, a);
                    newSmem.m = j;

                    if((newSmem.s < min_intv_array[i]) && ((smem.n - smem.m + 1) >= minSeedLen))
                    {
                        cur_j = j;

                        matchArray[numTotalSmem++] = smem;
                        break;
                    }
                    if((newSmem.s >= min_intv_array[i]) && (newSmem.s != curr_s))
                    {
                        curr_s = newSmem.s;
                        prev[numCurr++] = newSmem;
#ifdef ENABLE_PREFETCH
                        _mm_prefetch((const char *)(&cp_occ[(newSmem.k) >> CP_SHIFT]), _MM_HINT_T0);
                        _mm_prefetch((const char *)(&cp_occ[(newSmem.k + newSmem.s) >> CP_SHIFT]), _MM_HINT_T0);
#endif
                        break;
                    }
                }
                p++;
                for(; p < numPrev; p++)
                {
                    SMEM smem = prev[p];

                    SMEM newSmem = backwardExt(smem, a);
                    newSmem.m = j;


                    if((newSmem.s >= min_intv_array[i]) && (newSmem.s != curr_s))
                    {
                        curr_s = newSmem.s;
                        prev[numCurr++] = newSmem;
#ifdef ENABLE_PREFETCH
                        _mm_prefetch((const char *)(&cp_occ[(newSmem.k) >> CP_SHIFT]), _MM_HINT_T0);
                        _mm_prefetch((const char *)(&cp_occ[(newSmem.k + newSmem.s) >> CP_SHIFT]), _MM_HINT_T0);
#endif
                    }
                }
                numPrev = numCurr;
                if(numCurr == 0)
                {
                    break;
                }
            }
            if(numPrev != 0)
            {
                SMEM smem = prev[0];
                if(((smem.n - smem.m + 1) >= minSeedLen))
                {

                    matchArray[numTotalSmem++] = smem;
                }
                numPrev = 0;
            }
        }
        query_pos_array[i] = next_x;
    }
    (*__numTotalSmem) = numTotalSmem;
}

void FMI_search::getSMEMsAllPosOneThread(uint8_t *enc_qdb,
                                         int32_t *min_intv_array,
                                         int32_t *rid_array,
                                         int32_t numReads,
                                         int32_t batch_size,
                                         const bseq1_t *seq_,
                                         int32_t *query_cum_len_ar,
                                         int32_t max_readlength,
                                         int32_t minSeedLen,
                                         SMEM *matchArray,
                                         int64_t *__numTotalSmem)
{
    int16_t *query_pos_array = (int16_t *)_mm_malloc(numReads * sizeof(int16_t), 64);
    
    int32_t i;
    for(i = 0; i < numReads; i++)
        query_pos_array[i] = 0;

    int32_t numActive = numReads;
    (*__numTotalSmem) = 0;

    do
    {
        int32_t head = 0;
        int32_t tail = 0;
        for(head = 0; head < numActive; head++)
        {
            int readlength = seq_[rid_array[head]].l_seq;
            if(query_pos_array[head] < readlength)
            {
                rid_array[tail] = rid_array[head];
                query_pos_array[tail] = query_pos_array[head];
                min_intv_array[tail] = min_intv_array[head];
                tail++;             
            }               
        }
        getSMEMsOnePosOneThread(enc_qdb,
                                query_pos_array,
                                min_intv_array,
                                rid_array,
                                tail,
                                batch_size,
                                seq_,
                                query_cum_len_ar,
                                max_readlength,
                                minSeedLen,
                                matchArray,
                                __numTotalSmem);
        numActive = tail;
    } while(numActive > 0);

    _mm_free(query_pos_array);
}

int64_t FMI_search::bwtSeedStrategyAllPosOneThread(uint8_t *enc_qdb,
                                                   int32_t *max_intv_array,
                                                   int32_t numReads,
                                                   const bseq1_t *seq_,
                                                   int32_t *query_cum_len_ar,
                                                   int32_t minSeedLen,
                                                   SMEM *matchArray)
{
    int32_t i;

    int64_t numTotalSeed = 0;

    for(i = 0; i < numReads; i++)
    {
        int readlength = seq_[i].l_seq;
        int16_t x = 0;
        while(x < readlength)
        {
            int next_x = x + 1;

            // Forward search
            SMEM smem;
            smem.rid = i;
            smem.m = x;
            smem.n = x;
            
            int offset = query_cum_len_ar[i];
            uint8_t a = enc_qdb[offset + x];
            // uint8_t a = enc_qdb[i * readlength + x];

            if(a < 4)
            {
                smem.k = count[a];
                smem.l = count[3 - a];
                smem.s = count[a+1] - count[a];


                int j;
                for(j = x + 1; j < readlength; j++)
                {
                    next_x = j + 1;
                    // a = enc_qdb[i * readlength + j];
                    a = enc_qdb[offset + j];
                    if(a < 4)
                    {
                        SMEM smem_ = smem;

                        // Forward extension is backward extension with the BWT of reverse complement
                        smem_.k = smem.l;
                        smem_.l = smem.k;
                        SMEM newSmem_ = backwardExt(smem_, 3 - a);
                        //SMEM smem = backwardExt(smem, 3 - a);
                        //smem.n = j;
                        SMEM newSmem = newSmem_;
                        newSmem.k = newSmem_.l;
                        newSmem.l = newSmem_.k;
                        newSmem.n = j;
                        smem = newSmem;
#ifdef ENABLE_PREFETCH
                        _mm_prefetch((const char *)(&cp_occ[(smem.k) >> CP_SHIFT]), _MM_HINT_T0);
                        _mm_prefetch((const char *)(&cp_occ[(smem.l) >> CP_SHIFT]), _MM_HINT_T0);
#endif


                        if((smem.s < max_intv_array[i]) && ((smem.n - smem.m + 1) >= minSeedLen))
                        {

                            if(smem.s > 0)
                            {
                                matchArray[numTotalSeed++] = smem;
                            }
                            break;
                        }
                    }
                    else
                    {

                        break;
                    }
                }

            }
            x = next_x;
        }
    }
    return numTotalSeed;
}


void FMI_search::getSMEMs(uint8_t *enc_qdb,
        int32_t numReads,
        int32_t batch_size,
        int32_t readlength,
        int32_t minSeedLen,
        int32_t nthreads,
        SMEM *matchArray,
        int64_t *numTotalSmem)
{
    SMEM *prevArray = (SMEM *)_mm_malloc(nthreads * readlength * sizeof(SMEM), 64);
    SMEM *currArray = (SMEM *)_mm_malloc(nthreads * readlength * sizeof(SMEM), 64);


// #pragma omp parallel num_threads(nthreads)
    {
        int tid = 0; //omp_get_thread_num();   // removed omp
        numTotalSmem[tid] = 0;
        SMEM *myPrevArray = prevArray + tid * readlength;
        SMEM *myCurrArray = prevArray + tid * readlength;

        int32_t perThreadQuota = (numReads + (nthreads - 1)) / nthreads;
        int32_t first = tid * perThreadQuota;
        int32_t last  = (tid + 1) * perThreadQuota;
        if(last > numReads) last = numReads;
        SMEM *myMatchArray = matchArray + first * readlength;

        uint32_t i;
        // Perform SMEM for original reads
        for(i = first; i < last; i++)
        {
            int x = readlength - 1;
            int numPrev = 0;
            int numSmem = 0;

            while (x >= 0)
            {
                // Forward search
                SMEM smem;
                smem.rid = i;
                smem.m = x;
                smem.n = x;
                uint8_t a = enc_qdb[i * readlength + x];

                if(a > 3)
                {
                    x--;
                    continue;
                }
                smem.k = count[a];
                smem.l = count[3 - a];
                smem.s = count[a+1] - count[a];

                int j;
                for(j = x + 1; j < readlength; j++)
                {
                    a = enc_qdb[i * readlength + j];
                    if(a < 4)
                    {
                        SMEM smem_ = smem;

                        // Forward extension is backward extension with the BWT of reverse complement
                        smem_.k = smem.l;
                        smem_.l = smem.k;
                        SMEM newSmem_ = backwardExt(smem_, 3 - a);
                        SMEM newSmem = newSmem_;
                        newSmem.k = newSmem_.l;
                        newSmem.l = newSmem_.k;
                        newSmem.n = j;

                        if(newSmem.s != smem.s)
                        {
                            myPrevArray[numPrev] = smem;
                            numPrev++;
                        }
                        smem = newSmem;
                        if(newSmem.s == 0)
                        {
                            break;
                        }
                    }
                    else
                    {
                        myPrevArray[numPrev] = smem;
                        numPrev++;
                        break;
                    }
                }
                if(smem.s != 0)
                {
                    myPrevArray[numPrev++] = smem;
                }

                SMEM *curr, *prev;
                prev = myPrevArray;
                curr = myCurrArray;

                int p;
                for(p = 0; p < (numPrev/2); p++)
                {
                    SMEM temp = prev[p];
                    prev[p] = prev[numPrev - p - 1];
                    prev[numPrev - p - 1] = temp;
                }

                int next_x = x - 1;

                // Backward search
                int cur_j = readlength;
                for(j = x - 1; j >= 0; j--)
                {
                    int numCurr = 0;
                    int curr_s = -1;
                    a = enc_qdb[i * readlength + j];
                    //printf("a = %d\n", a);
                    if(a > 3)
                    {
                        next_x = j - 1;
                        break;
                    }
                    for(p = 0; p < numPrev; p++)
                    {
                        SMEM smem = prev[p];
                        SMEM newSmem = backwardExt(smem, a);
                        newSmem.m = j;

                        if(newSmem.s == 0)
                        {
                            if((numCurr == 0) && (j < cur_j))
                            {
                                cur_j = j;
                                if((smem.n - smem.m + 1) >= minSeedLen)
                                    myMatchArray[numTotalSmem[tid] + numSmem++] = smem;
                            }
                        }
                        if((newSmem.s != 0) && (newSmem.s != curr_s))
                        {
                            curr_s = newSmem.s;
                            curr[numCurr++] = newSmem;
                        }
                    }
                    SMEM *temp = prev;
                    prev = curr;
                    curr = temp;
                    numPrev = numCurr;
                    if(numCurr == 0)
                    {
                        next_x = j;
                        break;
                    }
                    else
                    {
                        next_x = j - 1;
                    }
                }
                if(numPrev != 0)
                {
                    SMEM smem = prev[0];
                    if((smem.n - smem.m + 1) >= minSeedLen)
                        myMatchArray[numTotalSmem[tid] + numSmem++] = smem;
                    numPrev = 0;
                }
                x = next_x;
            }
            numTotalSmem[tid] += numSmem;
        }
    }

    _mm_free(prevArray);
    _mm_free(currArray);
}


int compare_smem(const void *a, const void *b)
{
    SMEM *pa = (SMEM *)a;
    SMEM *pb = (SMEM *)b;

    if(pa->rid < pb->rid)
        return -1;
    if(pa->rid > pb->rid)
        return 1;

    if(pa->m < pb->m)
        return -1;
    if(pa->m > pb->m)
        return 1;
    if(pa->n > pb->n)
        return -1;
    if(pa->n < pb->n)
        return 1;
    return 0;
}

void FMI_search::sortSMEMs(SMEM *matchArray,
        int64_t numTotalSmem[],
        int32_t numReads,
        int32_t readlength,
        int nthreads)
{
    int tid;
    int32_t perThreadQuota = (numReads + (nthreads - 1)) / nthreads;
    for(tid = 0; tid < nthreads; tid++)
    {
        int32_t first = tid * perThreadQuota;
        SMEM *myMatchArray = matchArray + first * readlength;
        qsort(myMatchArray, numTotalSmem[tid], sizeof(SMEM), compare_smem);
    }
}


SMEM FMI_search::backwardExt(SMEM smem, uint8_t a)
{
    //beCalls++;
    uint8_t b;

    int64_t k[4], l[4], s[4];
    for(b = 0; b < 4; b++)
    {
        int64_t sp = (int64_t)(smem.k);
        int64_t ep = (int64_t)(smem.k) + (int64_t)(smem.s);
        GET_OCC_UNCOMPRESSED(sp, b, occ_sp);
        GET_OCC_UNCOMPRESSED(ep, b, occ_ep);
        //printf("GET OCC input: sp = %ld, b = %ld, ep = %ld; output: occ_sp = %ld, occ_ep = %ld \n", sp, b, ep, occ_sp, occ_ep);
        k[b] = count[b] + occ_sp;
        s[b] = occ_ep - occ_sp;
    }

    int64_t sentinel_offset = 0;
    if((smem.k <= sentinel_index) && ((smem.k + smem.s) > sentinel_index)) sentinel_offset = 1;
    l[3] = smem.l + sentinel_offset;
    l[2] = l[3] + s[3];
    l[1] = l[2] + s[2];
    l[0] = l[1] + s[1];

    smem.k = k[a];
    smem.l = l[a];
    smem.s = s[a];
    return smem;
}

int64_t FMI_search::get_sa_entry(int64_t pos)
{
    int64_t sa_entry = sa_ms_byte[pos];
    sa_entry = sa_entry << 32;
    sa_entry = sa_entry + sa_ls_word[pos];
    return sa_entry;
}

void FMI_search::get_sa_entries(int64_t *posArray, int64_t *coordArray, uint32_t count, int32_t nthreads)
{
    uint32_t i;
// #pragma omp parallel for num_threads(nthreads)
    for(i = 0; i < count; i++)
    {
        int64_t pos = posArray[i];
        int64_t sa_entry = sa_ms_byte[pos];
        sa_entry = sa_entry << 32;
        sa_entry = sa_entry + sa_ls_word[pos];
        //_mm_prefetch((const char *)(sa_ms_byte + pos + SAL_PFD), _MM_HINT_T0);
        coordArray[i] = sa_entry;
    }
}

void FMI_search::get_sa_entries(SMEM *smemArray, int64_t *coordArray, int32_t *coordCountArray, uint32_t count, int32_t max_occ)
{
    uint32_t i;
    int32_t totalCoordCount = 0;
    for(i = 0; i < count; i++)
    {
        int32_t c = 0;
        SMEM smem = smemArray[i];
        int64_t hi = smem.k + smem.s;
        int64_t step = (smem.s > max_occ) ? smem.s / max_occ : 1;
        int64_t j;
        for(j = smem.k; (j < hi) && (c < max_occ); j+=step, c++)
        {
            int64_t pos = j;
            int64_t sa_entry = sa_ms_byte[pos];
            sa_entry = sa_entry << 32;
            sa_entry = sa_entry + sa_ls_word[pos];
            //_mm_prefetch((const char *)(sa_ms_byte + pos + SAL_PFD * step), _MM_HINT_T0);
            coordArray[totalCoordCount + c] = sa_entry;
        }
        coordCountArray[i] = c;
        totalCoordCount += c;
    }
}

// sa_compression
int64_t FMI_search::get_sa_entry_compressed(int64_t pos, int tid)
{
    if ((pos & SA_COMPX_MASK) == 0) { 
        // if pos is a multiple of compression factor, just retrieve the entry from memory
        #if  SA_COMPRESSION
        int64_t sa_entry = sa_ms_byte[pos >> SA_COMPX];  // this looks like the most significant byte
        #else
        int64_t sa_entry = sa_ms_byte[pos];     // simulation
        #endif
        
        sa_entry = sa_entry << 32;
        
        #if  SA_COMPRESSION
        sa_entry = sa_entry + sa_ls_word[pos >> SA_COMPX];  // least significant word (32 bit)
        #else
        sa_entry = sa_entry + sa_ls_word[pos];   // simulation
        #endif
        
        return sa_entry;        
    }
    else {
        // else, the SA entry is not stored on memory and must be computed

        // tprof[MEM_CHAIN][tid] ++;
        int64_t offset = 0; 
        int64_t sp = pos;
        while(true)
        {
            int64_t occ_id_pp_ = sp >> CP_SHIFT;
            int64_t y_pp_ = CP_BLOCK_SIZE - (sp & CP_MASK) - 1; 
            uint8_t bwt_char = cp_occ[occ_id_pp_].bwt_char;

            if (bwt_char == 4) {
                return offset;
            }

            GET_OCC_UNCOMPRESSED(sp, bwt_char, occ_sp);
            // Occ[sp, BWT[sp]] , or, the rank of BWT[sp]

            // Go to the BW matrix row that start with BWT[sp] of that rank
            sp = count[bwt_char] + occ_sp;
            
            offset ++;
            // tprof[ALIGN1][tid] ++;
            // Stop until we reach a SA element row that is not compressed. Use offset to retrieve the request SA entry
            if ((sp & SA_COMPX_MASK) == 0) break;
        }
        // assert((reference_seq_len >> SA_COMPX) - 1 >= (sp >> SA_COMPX));
        #if  SA_COMPRESSION
        int64_t sa_entry = sa_ms_byte[sp >> SA_COMPX];
        #else
        int64_t sa_entry = sa_ms_byte[sp];      // simultion
        #endif
        
        sa_entry = sa_entry << 32;

        #if  SA_COMPRESSION
        sa_entry = sa_entry + sa_ls_word[sp >> SA_COMPX];
        #else
        sa_entry = sa_entry + sa_ls_word[sp];      // simulation
        #endif
        
        sa_entry += offset;
        return sa_entry;
    }
}

void FMI_search::get_sa_entries(SMEM *smemArray, int64_t *coordArray, int32_t *coordCountArray, uint32_t count, int32_t max_occ, int tid)
{
    
    uint32_t i;
    int32_t totalCoordCount = 0;
    for(i = 0; i < count; i++)
    {
        int32_t c = 0;
        SMEM smem = smemArray[i];
        int64_t hi = smem.k + smem.s;
        int64_t step = (smem.s > max_occ) ? smem.s / max_occ : 1;
        int64_t j;
        for(j = smem.k; (j < hi) && (c < max_occ); j+=step, c++)
        {
            int64_t pos = j;
            int64_t sa_entry = get_sa_entry_compressed(pos, tid);
            coordArray[totalCoordCount + c] = sa_entry;
        }
        // coordCountArray[i] = c;
        *coordCountArray += c;
        totalCoordCount += c;
    }
}

// SA_COPMRESSION w/ PREFETCH
int64_t FMI_search::call_one_step(int64_t pos, int64_t &sa_entry, int64_t &offset)
{
    if ((pos & SA_COMPX_MASK) == 0) {        
        sa_entry = sa_ms_byte[pos >> SA_COMPX];        
        sa_entry = sa_entry << 32;        
        sa_entry = sa_entry + sa_ls_word[pos >> SA_COMPX];        
        // return sa_entry;
        return 1;
    }
    else {
        // int64_t offset = 0; 
        int64_t sp = pos;

        int64_t occ_id_pp_ = sp >> CP_SHIFT;
        int64_t y_pp_ = CP_BLOCK_SIZE - (sp & CP_MASK) - 1; 
        uint8_t bwt_char = cp_occ[occ_id_pp_].bwt_char;
        uint8_t b;

        if (bwt_char == 4) {
            sa_entry = 0;
            return 1;
        }
        
        GET_OCC_UNCOMPRESSED(sp, bwt_char, occ_sp);
        
        sp = count[b] + occ_sp;
        
        offset ++;
        if ((sp & SA_COMPX_MASK) == 0) {
    
            sa_entry = sa_ms_byte[sp >> SA_COMPX];        
            sa_entry = sa_entry << 32;
            sa_entry = sa_entry + sa_ls_word[sp >> SA_COMPX];
            
            sa_entry += offset;
            // return sa_entry;
            return 1;
        }
        else {
            sa_entry = sp;
            return 0;
        }
    } // else
}

void FMI_search::get_sa_entries_prefetch(SMEM *smemArray, int64_t *coordArray,
                                         int64_t *coordCountArray, int64_t count,
                                         const int32_t max_occ, int tid, int64_t &id_)
{
    
    // uint32_t i;
    int32_t totalCoordCount = 0;
    int32_t mem_lim = 0, id = 0;
    
    for(int i = 0; i < count; i++)
    {
        int32_t c = 0;
        SMEM smem = smemArray[i];
        mem_lim += smem.s;
    }

    int64_t *pos_ar = (int64_t *) _mm_malloc( mem_lim * sizeof(int64_t), 64);
    int64_t *map_ar = (int64_t *) _mm_malloc( mem_lim * sizeof(int64_t), 64);

    for(int i = 0; i < count; i++)
    {
        int32_t c = 0;
        SMEM smem = smemArray[i];
        int64_t hi = smem.k + smem.s;
        int64_t step = (smem.s > max_occ) ? smem.s / max_occ : 1;
        int64_t j;
        for(j = smem.k; (j < hi) && (c < max_occ); j+=step, c++)
        {
            int64_t pos = j;
             pos_ar[id]  = pos;
             map_ar[id++] = totalCoordCount + c;
            // int64_t sa_entry = get_sa_entry_compressed(pos, tid);
            // coordArray[totalCoordCount + c] = sa_entry;
        }
        //coordCountArray[i] = c;
        *coordCountArray += c;
        totalCoordCount += c;
    }
    
    id_ += id;
    
    const int32_t sa_batch_size = 20;
    int64_t working_set[sa_batch_size], map_pos[sa_batch_size];;
    int64_t offset[sa_batch_size] = {-1};
    
    int i = 0, j = 0;    
    while(i<id && j<sa_batch_size)
    {
        int64_t pos =  pos_ar[i];
        working_set[j] = pos;
        map_pos[j] = map_ar[i];
        offset[j] = 0;
        
        if (pos & SA_COMPX_MASK == 0) {
            _mm_prefetch(&sa_ms_byte[pos >> SA_COMPX], _MM_HINT_T0);
            _mm_prefetch(&sa_ls_word[pos >> SA_COMPX], _MM_HINT_T0);
        }
        else {
            int64_t occ_id_pp_ = pos >> CP_SHIFT;
            _mm_prefetch(&cp_occ[occ_id_pp_], _MM_HINT_T0);
        }
        i++;
        j++;
    }
        
    int lim = j, all_quit = 0;
    while (all_quit < id)
    {
        
        for (int k=0; k<lim; k++)
        {
            int64_t sp = 0, pos = 0;
            bool quit;
            if (offset[k] >= 0) {
                quit = call_one_step(working_set[k], sp, offset[k]);
            }
            else
                continue;
            
            if (quit) {
                coordArray[map_pos[k]] = sp;
                all_quit ++;
                
                if (i < id)
                {
                    pos = pos_ar[i];
                    working_set[k] = pos;
                    map_pos[k] = map_ar[i++];
                    offset[k] = 0;
                    
                    if (pos & SA_COMPX_MASK == 0) {
                        _mm_prefetch(&sa_ms_byte[pos >> SA_COMPX], _MM_HINT_T0);
                        _mm_prefetch(&sa_ls_word[pos >> SA_COMPX], _MM_HINT_T0);
                    }
                    else {
                        int64_t occ_id_pp_ = pos >> CP_SHIFT;
                        _mm_prefetch(&cp_occ[occ_id_pp_], _MM_HINT_T0);
                    }
                }
                else
                    offset[k] = -1;
            }
            else {
                working_set[k] = sp;
                if (sp & SA_COMPX_MASK == 0) {
                    _mm_prefetch(&sa_ms_byte[sp >> SA_COMPX], _MM_HINT_T0);
                    _mm_prefetch(&sa_ls_word[sp >> SA_COMPX], _MM_HINT_T0);
                }
                else {
                    int64_t occ_id_pp_ = sp >> CP_SHIFT;
                    _mm_prefetch(&cp_occ[occ_id_pp_], _MM_HINT_T0);
                }                
            }
        }
    }
    
    _mm_free(pos_ar);
    _mm_free(map_ar);
}
