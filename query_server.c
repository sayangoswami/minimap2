//
// Created by sayan on 25.04.24.
//

#include "minimap.h"
#include "mmpriv.h"
#include "ketopt.h"
#include "kvec.h"
#include "kalloc.h"
#include "sdust.h"
#include "bseq.h"
#include "khash.h"
#include <errno.h>
#include <signal.h>

/** generic macros */
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define alignup(n, a) (((n) + (a)-1) & ~((a)-1))

/** string macros */
#define streq(s1, s2) (!strcmp((s1), (s2)))
#define str_endswith(str, suffix) (streq((suffix), (str) + (strlen(str) - strlen(suffix))))

#define ERR(fmt, ...) ({ \
    fprintf(stderr, "ERROR: " fmt " at %s:%i\n",	\
	    ##__VA_ARGS__, __FILE__, __LINE__);	\
    raise(SIGTRAP);									\
})

#define expect(expression) if (!(expression)) ERR("Expected " #expression "")

typedef struct {
    int n_threads, n_parts, old_best_n;
    char *fnw, *rg, *junc_bed, *s, *alt_list;
} mm_genopt_t;

typedef struct {
    mm_mapopt_t opt;
    mm_idxopt_t ipt;
    mm_genopt_t gpt;
    mm_idx_t *idx;
    mm_tbuf_t *b;
} sdata_t;

typedef struct {
    int num_seeds, max_count, num_chains;
    uint64_t num_anchors, sum_scores,
            sum_len_seeded_match,
            sum_len_alignment_blk,
            sum_best_alt_map_score,
            sum_num_suboptimal_mappings,
            sum_dp_max, sum_dp_max2,
            sum_dp_scores, sum_initial_chaining_scores;
    double sum_est_err;
} result_t;

extern void collect_minimizers(void *km, const mm_mapopt_t *opt, const mm_idx_t *mi, int n_segs,
                               const int *qlens, const char **seqs, mm128_v *mv);
extern mm128_t *collect_seed_hits_heap(void *km, const mm_mapopt_t *opt, int max_occ, const mm_idx_t *mi,
                                       const char *qname, const mm128_v *mv, int qlen, int64_t *n_a, int *rep_len,
                                       int *n_mini_pos, uint64_t **mini_pos);
extern mm128_t *collect_seed_hits(void *km, const mm_mapopt_t *opt, int max_occ, const mm_idx_t *mi,
                           const char *qname, const mm128_v *mv, int qlen, int64_t *n_a, int *rep_len,
                           int *n_mini_pos, uint64_t **mini_pos);
extern void chain_post(const mm_mapopt_t *opt, int max_chain_gap_ref, const mm_idx_t *mi,
                       void *km, int qlen, int n_segs, const int *qlens, int *n_regs, mm_reg1_t *regs, mm128_t *a);
extern mm_reg1_t *align_regs(const mm_mapopt_t *opt, const mm_idx_t *mi, void *km, int qlen,
                             const char *seq, int *n_regs, mm_reg1_t *regs, mm128_t *a);


static void* mmidx_build2(char *ref, const int n_threads, const uint64_t batch_size_gb) {
    mm_verbose = 3;
    sdata_t *sdata = malloc(sizeof(sdata_t));
    sdata->gpt.n_threads = n_threads;
    sdata->gpt.old_best_n = -1;
    sdata->gpt.fnw = 0, sdata->gpt.rg = 0, sdata->gpt.junc_bed = 0, sdata->gpt.alt_list = 0;
    int rc = mm_set_opt(0, &sdata->ipt, &sdata->opt);
    expect(rc == 0);
    sdata->ipt.batch_size = MAX(sdata->ipt.batch_size, batch_size_gb<<30);
    mm_idx_reader_t *idx_rdr = mm_idx_reader_open(ref, &sdata->ipt, NULL);
    if (!idx_rdr) ERR("failed to open file %s (%s)", ref, strerror(errno));

    if ((sdata->idx = mm_idx_reader_read(idx_rdr, sdata->gpt.n_threads)) != 0) {
        int ret;
        if ((sdata->opt.flag & MM_F_CIGAR) && (sdata->idx->flag & MM_I_NO_SEQ)) {
            fprintf(stderr, "[ERROR] the prebuilt index doesn't contain sequences.\n");
            mm_idx_destroy(sdata->idx);
            mm_idx_reader_close(idx_rdr);
            return NULL;
        }
        if ((sdata->opt.flag & MM_F_OUT_SAM) && idx_rdr->n_parts == 1) {
            ERR("Not Implemented");
        }
        if (mm_verbose >= 3)
            fprintf(stderr, "[M::%s::%.3f*%.2f] loaded/built the index for %d target sequence(s)\n",
                    __func__, realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0), sdata->idx->n_seq);
        mm_mapopt_update(&sdata->opt, sdata->idx);
        if (mm_verbose >= 3) mm_idx_stat(sdata->idx);
        if (sdata->gpt.junc_bed) mm_idx_bed_read(sdata->idx, sdata->gpt.junc_bed, 1);
        if (sdata->gpt.alt_list) mm_idx_alt_read(sdata->idx, sdata->gpt.alt_list);

        sdata->b = mm_tbuf_init();
    }
    else ERR("Unexpectedly reached end-of-file.");

    return sdata;
}

static int mm_map_crude2(const mm_idx_t *mi, const int qlen, const char *seq, mm_tbuf_t *b,
                         const mm_mapopt_t *opt, const char *qname, result_t *result)
{
    int i, j, rep_len, qlen_sum = qlen, n_regs0, n_mini_pos;
    int max_chain_gap_qry, max_chain_gap_ref, is_splice = !!(opt->flag & MM_F_SPLICE), is_sr = !!(opt->flag & MM_F_SR);
    uint32_t hash;
    int64_t n_a;
    uint64_t *u, *mini_pos;
    mm128_t *a;
    mm128_v mv = {0,0,0};
    mm_reg1_t *regs0;
    km_stat_t kmst;
    float chn_pen_gap, chn_pen_skip;

    if (qlen_sum == 0) return -1;
    if (opt->max_qlen > 0 && qlen_sum > opt->max_qlen) return -1;

    hash  = qname && !(opt->flag & MM_F_NO_HASH_NAME)? __ac_X31_hash_string(qname) : 0;
    hash ^= __ac_Wang_hash(qlen_sum) + __ac_Wang_hash(opt->seed);
    hash  = __ac_Wang_hash(hash);

    collect_minimizers(b->km, opt, mi, 1, &qlen, &seq, &mv);
    if (opt->q_occ_frac > 0.0f) mm_seed_mz_flt(b->km, &mv, opt->mid_occ, opt->q_occ_frac);
    if (opt->flag & MM_F_HEAP_SORT) a = collect_seed_hits_heap(b->km, opt, opt->mid_occ, mi, qname, &mv, qlen_sum, &n_a, &rep_len, &n_mini_pos, &mini_pos);
    else a = collect_seed_hits(b->km, opt, opt->mid_occ, mi, qname, &mv, qlen_sum, &n_a, &rep_len, &n_mini_pos, &mini_pos);

    /********************** STAGE 1 **********************/
    if (LATENCY_STAGE == 1) {
        result->num_seeds = n_a;
        int max_count = 0, count = 0;
        for (i = 0; i < n_a-1; ++i) {
            if (a[i].x >> 32 == a[i+1].x >> 32) count++;
            else if (count > max_count) max_count = count, count = 0;
        }
        result->max_count = max_count;
        kfree(b->km, mv.a);
        kfree(b->km, a);
        kfree(b->km, mini_pos);
        return 0;
    }
    /****************************************************/

    // set max chaining gap on the query and the reference sequence
    if (is_sr)
        max_chain_gap_qry = qlen_sum > opt->max_gap? qlen_sum : opt->max_gap;
    else max_chain_gap_qry = opt->max_gap;
    if (opt->max_gap_ref > 0) {
        max_chain_gap_ref = opt->max_gap_ref; // always honor mm_mapopt_t::max_gap_ref if set
    } else if (opt->max_frag_len > 0) {
        max_chain_gap_ref = opt->max_frag_len - qlen_sum;
        if (max_chain_gap_ref < opt->max_gap) max_chain_gap_ref = opt->max_gap;
    } else max_chain_gap_ref = opt->max_gap;

    chn_pen_gap  = opt->chain_gap_scale * 0.01 * mi->k;
    chn_pen_skip = opt->chain_skip_scale * 0.01 * mi->k;
    if (opt->flag & MM_F_RMQ) {
        a = mg_lchain_rmq(opt->max_gap, opt->rmq_inner_dist, opt->bw, opt->max_chain_skip, opt->rmq_size_cap, opt->min_cnt, opt->min_chain_score,
                          chn_pen_gap, chn_pen_skip, n_a, a, &n_regs0, &u, b->km);
    } else {
        a = mg_lchain_dp(max_chain_gap_ref, max_chain_gap_qry, opt->bw, opt->max_chain_skip, opt->max_chain_iter, opt->min_cnt, opt->min_chain_score,
                         chn_pen_gap, chn_pen_skip, is_splice, 1, n_a, a, &n_regs0, &u, b->km);
    }

    /********************** STAGE 2 **********************/
    if (LATENCY_STAGE == 2) {
        result->num_chains = n_regs0;
        result->num_anchors = result->sum_scores = result->max_count = 0;
        for (int i = 0; i < n_regs0; ++i) {
            result->sum_scores += ((u[i] >> 32) & 0xffffffff), result->num_anchors += (u[i] & 0xffffffff);
        }
        int count = 0;
        for (int i = 0; i < result->num_anchors-1; ++i) {
            if (a[i].x >> 32 == a[i+1].x >> 32) count++;
            else if (count > result->max_count) result->max_count = count, count = 0;
        }

        goto CLEANUP; //< number of chains
    }
    /****************************************************/

    if (opt->bw_long > opt->bw && (opt->flag & (MM_F_SPLICE|MM_F_SR|MM_F_NO_LJOIN)) == 0 && n_regs0 > 1) { // re-chain/long-join for long sequences
        int32_t st = (int32_t)a[0].y, en = (int32_t)a[(int32_t)u[0] - 1].y;
        if (qlen_sum - (en - st) > opt->rmq_rescue_size || en - st > qlen_sum * opt->rmq_rescue_ratio) {
            int32_t i;
            for (i = 0, n_a = 0; i < n_regs0; ++i) n_a += (int32_t)u[i];
            kfree(b->km, u);
            radix_sort_128x(a, a + n_a);
            a = mg_lchain_rmq(opt->max_gap, opt->rmq_inner_dist, opt->bw_long, opt->max_chain_skip, opt->rmq_size_cap, opt->min_cnt, opt->min_chain_score,
                              chn_pen_gap, chn_pen_skip, n_a, a, &n_regs0, &u, b->km);
        }
    } else if (opt->max_occ > opt->mid_occ && rep_len > 0 && !(opt->flag & MM_F_RMQ)) { // re-chain, mostly for short reads
        int rechain = 0;
        if (n_regs0 > 0) { // test if the best chain has all the segments
            int n_chained_segs = 1, max = 0, max_i = -1, max_off = -1, off = 0;
            for (i = 0; i < n_regs0; ++i) { // find the best chain
                if (max < (int)(u[i]>>32)) max = u[i]>>32, max_i = i, max_off = off;
                off += (uint32_t)u[i];
            }
            for (i = 1; i < (int32_t)u[max_i]; ++i) // count the number of segments in the best chain
                if ((a[max_off+i].y&MM_SEED_SEG_MASK) != (a[max_off+i-1].y&MM_SEED_SEG_MASK))
                    ++n_chained_segs;
            if (n_chained_segs < 1)
                rechain = 1;
        } else rechain = 1;
        if (rechain) { // redo chaining with a higher max_occ threshold
            kfree(b->km, a);
            kfree(b->km, u);
            kfree(b->km, mini_pos);
            if (opt->flag & MM_F_HEAP_SORT) a = collect_seed_hits_heap(b->km, opt, opt->max_occ, mi, qname, &mv, qlen_sum, &n_a, &rep_len, &n_mini_pos, &mini_pos);
            else a = collect_seed_hits(b->km, opt, opt->max_occ, mi, qname, &mv, qlen_sum, &n_a, &rep_len, &n_mini_pos, &mini_pos);
            a = mg_lchain_dp(max_chain_gap_ref, max_chain_gap_qry, opt->bw, opt->max_chain_skip, opt->max_chain_iter, opt->min_cnt, opt->min_chain_score,
                             chn_pen_gap, chn_pen_skip, is_splice, 1, n_a, a, &n_regs0, &u, b->km);
        }
    }
    b->frag_gap = max_chain_gap_ref;
    b->rep_len = rep_len;

    /********************** STAGE 3 **********************/
    if (LATENCY_STAGE == 3) {
        result->num_seeds = n_a;
        result->num_chains = n_regs0;
        result->num_anchors = result->sum_scores = 0;
        for (int i = 0; i < n_regs0; ++i)
            result->sum_scores += ((u[i]>>32) & 0xffffffff), result->num_anchors += (u[i] & 0xffffffff);
        goto CLEANUP; //< refined number of chains
    }
    /****************************************************/

    regs0 = mm_gen_regs(b->km, hash, qlen_sum, n_regs0, u, a, !!(opt->flag&MM_F_QSTRAND));

    if (mi->n_alt) {
        mm_mark_alt(mi, n_regs0, regs0);
        mm_hit_sort(b->km, &n_regs0, regs0, opt->alt_drop); // this step can be merged into mm_gen_regs(); will do if this shows up in profile
    }

    /********************** STAGE 4 **********************/
    if (LATENCY_STAGE == 4) {
        result->num_seeds = n_a;
        result->num_chains = n_regs0;
        result->num_anchors = result->sum_scores = 0;
        result->sum_len_seeded_match = result->sum_len_alignment_blk = 0;
        for (int i = 0; i < n_regs0; ++i) {
            result->sum_scores += ((u[i]>>32) & 0xffffffff), result->num_anchors += (u[i] & 0xffffffff);
            result->sum_len_seeded_match += regs0[i].mlen, result->sum_len_alignment_blk += regs0[i].blen;
        }
        goto CLEANUP; //< length of chains
    }
    /****************************************************/

    chain_post(opt, max_chain_gap_ref, mi, b->km, qlen_sum, 1, &qlen, &n_regs0, regs0, a);

    /********************** STAGE 5 **********************/
    if (LATENCY_STAGE == 5) {
        result->num_seeds = n_a;
        result->num_chains = n_regs0;
        result->num_anchors = result->sum_scores = 0;
        result->sum_len_seeded_match = result->sum_len_alignment_blk = 0;
        result->sum_best_alt_map_score = result->sum_num_suboptimal_mappings = result->sum_dp_max = result->sum_dp_max2 = 0;
        for (int i = 0; i < n_regs0; ++i) {
            result->sum_scores += ((u[i]>>32) & 0xffffffff), result->num_anchors += (u[i] & 0xffffffff);
            result->sum_len_seeded_match += regs0[i].mlen, result->sum_len_alignment_blk += regs0[i].blen;
            result->sum_num_suboptimal_mappings += regs0[i].n_sub;
            result->sum_best_alt_map_score += regs0[i].subsc;
            if (regs0[i].p) result->sum_dp_max += regs0[i].p->dp_max, result->sum_dp_max2 += regs0[i].p->dp_max2;
        }
        goto CLEANUP;
    }
    /****************************************************/

    if (!is_sr && !(opt->flag&MM_F_QSTRAND)) {
        mm_est_err(mi, qlen_sum, n_regs0, regs0, a, n_mini_pos, mini_pos);
        n_regs0 = mm_filter_strand_retained(n_regs0, regs0);
    }

    regs0 = align_regs(opt, mi, b->km, qlen, seq, &n_regs0, regs0, a);
    regs0 = (mm_reg1_t*)realloc(regs0, sizeof(*regs0) * n_regs0);
    mm_set_mapq(b->km, n_regs0, regs0, opt->min_chain_score, opt->a, rep_len, is_sr);

    /********************** STAGE 6 **********************/
    result->num_seeds = n_a;
    result->num_chains = n_regs0;
    result->num_anchors = result->sum_scores = 0;
    result->sum_len_seeded_match = result->sum_len_alignment_blk = 0;
    result->sum_best_alt_map_score = result->sum_num_suboptimal_mappings = result->sum_dp_max = result->sum_dp_max2 = 0;
    result->sum_est_err = 0.0f;
    result->sum_dp_scores = result->sum_initial_chaining_scores = 0;
    for (int i = 0; i < n_regs0; ++i) {
        result->sum_scores += ((u[i]>>32) & 0xffffffff), result->num_anchors += (u[i] & 0xffffffff);
        result->sum_len_seeded_match += regs0[i].mlen, result->sum_len_alignment_blk += regs0[i].blen;
        result->sum_num_suboptimal_mappings += regs0[i].n_sub;
        result->sum_best_alt_map_score += regs0[i].subsc;
        if (regs0[i].p) result->sum_dp_max += regs0[i].p->dp_max, result->sum_dp_max2 += regs0[i].p->dp_max2;
        result->sum_est_err += regs0[i].div;
        result->sum_dp_scores += regs0[i].score, result->sum_initial_chaining_scores += regs0[i].score0;
    }
    /****************************************************/

    CLEANUP:
    kfree(b->km, mv.a);
    kfree(b->km, a);
    kfree(b->km, u);
    kfree(b->km, mini_pos);

    if (b->km) {
        km_stat(b->km, &kmst);
        assert(kmst.n_blocks == kmst.n_cores); // otherwise, there is a memory leak
        if (kmst.largest > 1U << 28 || (opt->cap_kalloc > 0 && kmst.capacity > opt->cap_kalloc)) {
            km_destroy(b->km);
            b->km = km_init();
        }
    }

    /// originally this was handled by the caller
    for (j = 0; j < n_regs0; j++)
        if (LATENCY_STAGE > 5 && regs0[j].p)
            free(regs0[j].p);
    if (n_regs0)
        free(regs0);
    return 0;
}


void* mmidx_build(char *ref) {
    return mmidx_build2(ref, 32, 4);
}

void mm_get_counts(void *index, const char *seq, int len, u_long *stats, void *tbuf) {
    sdata_t *sdata = index;
    result_t result;
    mm_map_crude2(sdata->idx, len, seq, tbuf, &sdata->opt, "q", &result);

    switch (LATENCY_STAGE) {
        case 1:
            stats[0] = result.num_seeds;
            stats[1] = result.max_count;
            break;
        case 2:
        case 3:
            stats[0] = result.num_chains;
            stats[1] = result.num_anchors;
            stats[2] = result.sum_scores;
            break;
        case 4:
            stats[0] = result.num_chains;
            stats[1] = result.sum_len_seeded_match;
            stats[2] = result.sum_len_alignment_blk;
            break;
        case 5:
        case 6:
            ERR("Not implemented");
    }
}