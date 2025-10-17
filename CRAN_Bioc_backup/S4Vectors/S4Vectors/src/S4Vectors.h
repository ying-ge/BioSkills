#include "../inst/include/S4Vectors_defines.h"
#include <string.h>

#define INIT_STATIC_SYMBOL(NAME) \
{ \
	if (NAME ## _symbol == NULL) \
		NAME ## _symbol = install(# NAME); \
}


/* safe_arithm.c */

void _reset_ovflow_flag();

int _get_ovflow_flag();

int _safe_int_add(
	int x,
	int y
);

int _safe_int_subtract(
	int x,
	int y
);

int _safe_int_mult(
	int x,
	int y
);

int _as_int(
	const char *val,
	int val_len
);

long long int _safe_llint_add(
	long long int x,
	long long int y
);

long long int _safe_llint_subtract(
	long long int x,
	long long int y
);

long long int _safe_llint_mult(
	long long int x,
	long long int y
);


/* sort_utils.c */

SEXP test_sort_ushort_array(
	SEXP x,
	SEXP desc
);

void _sort_int_array(
	int *x,
	size_t nelt,
	int desc
);

void _get_order_of_int_array(
	const int *x,
	int nelt,
	int desc,
	int *out,
	int out_shift
);

int _sort_ints(
	int *base,
	int base_len,
	const int *x,
	int desc,
	int use_radix,
	unsigned short int *rxbuf1,
	int *rxbuf2
);

void _pcompare_int_pairs(
	const int *a1,
	const int *b1,
	int nelt1,
	const int *a2,
	const int *b2,
	int nelt2,
	int *out,
	int out_len,
	int with_warning
);

int _int_pairs_are_sorted(
	const int *a,
	const int *b,
	int nelt,
	int desc,
	int strict
);

void _get_order_of_int_pairs(
	const int *a,
	const int *b,
	int nelt,
	int a_desc,
	int b_desc,
	int *out,
	int out_shift
);

int _sort_int_pairs(
	int *base,
	int base_len,
	const int *a,
	const int *b,
	int a_desc,
	int b_desc,
	int use_radix,
	unsigned short int *rxbuf1,
	int *rxbuf2
);

void _get_matches_of_ordered_int_pairs(
	const int *a1,
	const int *b1,
	const int *o1,
	int nelt1,
	const int *a2,
	const int *b2,
	const int *o2,
	int nelt2,
	int nomatch,
	int *out,
	int out_shift
);

int _int_quads_are_sorted(
	const int *a,
	const int *b,
	const int *c,
	const int *d,
	int nelt,
	int desc,
	int strict
);

void _get_order_of_int_quads(
	const int *a,
	const int *b,
	const int *c,
	const int *d,
	int nelt,
	int a_desc,
	int b_desc,
	int c_desc,
	int d_desc,
	int *out,
	int out_shift
);

int _sort_int_quads(
	int *base,
	int base_len,
	const int *a,
	const int *b,
	const int *c,
	const int *d,
	int a_desc,
	int b_desc,
	int c_desc,
	int d_desc,
	int use_radix,
	unsigned short int *rxbuf1,
	int *rxbuf2
);

void _get_matches_of_ordered_int_quads(
	const int *a1,
	const int *b1,
	const int *c1,
	const int *d1,
	const int *o1,
	int nelt1,
	const int *a2,
	const int *b2,
	const int *c2,
	const int *d2,
	const int *o2,
	int nelt2,
	int nomatch,
	int *out,
	int out_shift
);


/* hash_utils.c */

struct htab _new_htab(int n);

int _get_hbucket_val(
	const struct htab *htab,
	int bucket_idx
);

void _set_hbucket_val(
	struct htab *htab,
	int bucket_idx,
	int val
);


/* AEbufs.c */

SEXP AEbufs_use_malloc(SEXP x);

size_t _increase_buflength(size_t buflength);

size_t _IntAE_get_nelt(const IntAE *ae);

size_t _IntAE_set_nelt(
	IntAE *ae,
	size_t nelt
);

void _IntAE_set_val(
	const IntAE *ae,
	int val
);

void _IntAE_extend(
	IntAE *ae,
	size_t new_buflength
);

void _IntAE_insert_at(
	IntAE *ae,
	size_t at,
	int val
);

IntAE *_new_IntAE(
	size_t buflength,
	size_t nelt,
	int val
);

void _IntAE_append(
	IntAE *ae,
	const int *newvals,
	size_t nnewval
);

void _IntAE_delete_at(
	IntAE *ae,
	size_t at,
	size_t nelt
);

void _IntAE_shift(
	const IntAE *ae,
	size_t offset,
	int shift
);

void _IntAE_sum_and_shift(
	const IntAE *ae1,
	const IntAE *ae2,
	int shift
);

void _IntAE_qsort(
	const IntAE *ae,
	size_t offset,
	int desc
);

void _IntAE_uniq(
	IntAE *ae,
	size_t offset
);

SEXP _new_INTEGER_from_IntAE(const IntAE *ae);

SEXP _new_LOGICAL_from_IntAE(const IntAE *ae);

IntAE *_new_IntAE_from_INTEGER(SEXP x);

IntAE *_new_IntAE_from_CHARACTER(
	SEXP x,
	int keyshift
);

size_t _IntAEAE_get_nelt(const IntAEAE *aeae);

size_t _IntAEAE_set_nelt(
	IntAEAE *aeae,
	size_t nelt
);

void _IntAEAE_extend(
	IntAEAE *aeae,
	size_t new_buflength
);

void _IntAEAE_insert_at(
	IntAEAE *aeae,
	size_t at,
	IntAE *ae
);

IntAEAE *_new_IntAEAE(
	size_t buflength,
	size_t nelt
);

void _IntAEAE_pappend(
	const IntAEAE *aeae1,
	const IntAEAE *aeae2
);

void _IntAEAE_shift(
	const IntAEAE *aeae,
	int shift
);

void _IntAEAE_sum_and_shift(
	const IntAEAE *aeae1,
	const IntAEAE *aeae2,
	int shift
);

SEXP _new_LIST_from_IntAEAE(
	const IntAEAE *aeae,
	int mode
);

IntAEAE *_new_IntAEAE_from_LIST(SEXP x);

SEXP _IntAEAE_toEnvir(
	const IntAEAE *aeae,
	SEXP envir,
	int keyshift
);

size_t _IntPairAE_get_nelt(const IntPairAE *ae);

size_t _IntPairAE_set_nelt(
	IntPairAE *ae,
	size_t nelt
);

void _IntPairAE_extend(
	IntPairAE *ae,
	size_t new_buflength
);

void _IntPairAE_insert_at(
	IntPairAE *ae,
	size_t at,
	int a,
	int b
);

IntPairAE *_new_IntPairAE(
	size_t buflength,
	size_t nelt
);

size_t _IntPairAEAE_get_nelt(const IntPairAEAE *aeae);

size_t _IntPairAEAE_set_nelt(
	IntPairAEAE *aeae,
	size_t nelt
);

void _IntPairAEAE_extend(
	IntPairAEAE *aeae,
	size_t new_buflength
);

void _IntPairAEAE_insert_at(
	IntPairAEAE *aeae,
	size_t at,
	IntPairAE *ae
);

IntPairAEAE *_new_IntPairAEAE(
	size_t buflength,
	size_t nelt
);

size_t _LLongAE_get_nelt(const LLongAE *ae);

size_t _LLongAE_set_nelt(
	LLongAE *ae,
	size_t nelt
);

void _LLongAE_set_val(
	const LLongAE *ae,
	long long val
);

void _LLongAE_extend(
	LLongAE *ae,
	size_t new_buflength
);

void _LLongAE_insert_at(
	LLongAE *ae,
	size_t at,
	long long val
);

LLongAE *_new_LLongAE(
	size_t buflength,
	size_t nelt,
	long long val
);

size_t _LLongAEAE_get_nelt(const LLongAEAE *aeae);

size_t _LLongAEAE_set_nelt(
	LLongAEAE *aeae,
	size_t nelt
);

void _LLongAEAE_extend(
	LLongAEAE *aeae,
	size_t new_buflength
);

void _LLongAEAE_insert_at(
	LLongAEAE *aeae,
	size_t at,
	LLongAE *ae
);

LLongAEAE *_new_LLongAEAE(
	size_t buflength,
	size_t nelt
);

size_t _DoubleAE_get_nelt(const DoubleAE *ae);

size_t _DoubleAE_set_nelt(
	DoubleAE *ae,
	size_t nelt
);

void _DoubleAE_set_val(
	const DoubleAE *ae,
	double val
);

void _DoubleAE_extend(
	DoubleAE *ae,
	size_t new_buflength
);

void _DoubleAE_insert_at(
	DoubleAE *ae,
	size_t at,
	double val
);

DoubleAE *_new_DoubleAE(
	size_t buflength,
	size_t nelt,
	double val
);

void _DoubleAE_append(
	DoubleAE *ae,
	const double *newvals,
	size_t nnewval
);

void _DoubleAE_delete_at(
	DoubleAE *ae,
	size_t at,
	size_t nelt
);

SEXP _new_NUMERIC_from_DoubleAE(const DoubleAE *ae);

DoubleAE *_new_DoubleAE_from_NUMERIC(SEXP x);

size_t _CharAE_get_nelt(const CharAE *ae);

size_t _CharAE_set_nelt(
	CharAE *ae,
	size_t nelt
);

void _CharAE_extend(
	CharAE *ae,
	size_t new_buflength
);

void _CharAE_insert_at(
	CharAE *ae,
	size_t at,
	char c
);

CharAE *_new_CharAE(size_t buflength);

CharAE *_new_CharAE_from_string(const char *string);

void _CharAE_append_string(
	CharAE *ae,
	const char *string
);

void _CharAE_delete_at(
	CharAE *ae,
	size_t at,
	size_t nelt
);

SEXP _new_CHARSXP_from_CharAE(const CharAE *ae);

SEXP _new_RAW_from_CharAE(const CharAE *ae);

SEXP _new_LOGICAL_from_CharAE(const CharAE *ae);

size_t _CharAEAE_get_nelt(const CharAEAE *aeae);

size_t _CharAEAE_set_nelt(
	CharAEAE *aeae,
	size_t nelt
);

void _CharAEAE_extend(
	CharAEAE *aeae,
	size_t new_buflength
);

void _CharAEAE_insert_at(
	CharAEAE *aeae,
	size_t at,
	CharAE *ae
);

CharAEAE *_new_CharAEAE(
	size_t buflength,
	size_t nelt
);

void _CharAEAE_append_string(
	CharAEAE *aeae,
	const char *string
);

SEXP _new_CHARACTER_from_CharAEAE(const CharAEAE *aeae);

SEXP AEbufs_free();


/* SEXP_utils.c */

const char *_get_classname(SEXP x);


/* anyMissing.c */

SEXP anyMissing(SEXP x);


/* LLint_class.c */

int _is_LLint(SEXP x);

SEXP make_RAW_from_NA_LLINT();

int sscan_llint(
	const char *s,
	int maxparse,
	long long int *val,
	int parse_dec
);

R_xlen_t _get_LLint_length(SEXP x);

long long int *_get_LLint_dataptr(SEXP x);

SEXP _alloc_LLint(const char *classname, R_xlen_t length);

SEXP new_LLint_from_LOGICAL(SEXP x);

SEXP new_LLint_from_INTEGER(SEXP x);

SEXP new_LLint_from_NUMERIC(SEXP x);

SEXP new_LLint_from_CHARACTER(SEXP x);

SEXP new_LOGICAL_from_LLint(SEXP x);

SEXP new_INTEGER_from_LLint(SEXP x);

SEXP new_NUMERIC_from_LLint(SEXP x);

SEXP new_CHARACTER_from_LLint(SEXP x);


/* subsetting_utils.c */

long long int _copy_vector_block(
	SEXP dest,
	long long int dest_offset,
	SEXP src,
	long long int src_offset,
	long long int block_nelt
);

int _copy_vector_positions(
	SEXP dest,
	int dest_offset,
	SEXP src,
	const int *pos,
	int npos
);

int _copy_vector_ranges(
	SEXP dest,
	int dest_offset,
	SEXP src,
	const int *start,
	const int *width,
	int nranges
);

SEXP _subset_vector_OR_factor_by_positions(
	SEXP x,
	const int *pos,
	int npos
);

SEXP _subset_vector_OR_factor_by_ranges(
	SEXP x,
	const int *start,
	const int *width,
	int nranges
);

SEXP vector_OR_factor_extract_positions(
	SEXP x,
	SEXP pos
);

SEXP vector_OR_factor_extract_ranges(
	SEXP x,
	SEXP start,
	SEXP width
);


/* vector_utils.c */

int _vector_memcmp(
	SEXP x1,
	int x1_offset,
	SEXP x2,
	int x2_offset,
	int nelt
);

SEXP sapply_NROW(SEXP x);

SEXP _list_as_data_frame(
	SEXP x,
	int nrow
);


/* logical_utils.c */

SEXP logical_sum(
	SEXP x,
	SEXP na_rm
);

SEXP logical2_sum(
	SEXP x,
	SEXP na_rm
);


/* integer_utils.c */

SEXP to_list_of_ints(
	SEXP x,
	SEXP sep
);

SEXP Integer_any_missing_or_outside(SEXP x, SEXP lower, SEXP upper);

SEXP Integer_diff_with_0(SEXP x);

SEXP Integer_diff_with_last(SEXP x, SEXP last);

SEXP Integer_order(
	SEXP x,
	SEXP decreasing,
	SEXP use_radix
);

int _check_integer_pairs(
	SEXP a,
	SEXP b,
	const int **a_p,
	const int **b_p,
	const char *a_argname,
	const char *b_argname
);

SEXP Integer_pcompare2(
	SEXP a1,
	SEXP b1,
	SEXP a2,
	SEXP b2
);

SEXP Integer_sorted2(
	SEXP a,
	SEXP b,
	SEXP decreasing,
	SEXP strictly
);

SEXP Integer_order2(
	SEXP a,
	SEXP b,
	SEXP decreasing,
	SEXP use_radix
);

SEXP Integer_match2_quick(
	SEXP a1,
	SEXP b1,
	SEXP a2,
	SEXP b2,
	SEXP nomatch
);

SEXP Integer_selfmatch2_quick(
	SEXP a,
	SEXP b
);

SEXP Integer_match2_hash(
	SEXP a1,
	SEXP b1,
	SEXP a2,
	SEXP b2,
	SEXP nomatch
);

SEXP Integer_selfmatch2_hash(
	SEXP a,
	SEXP b
);

int _check_integer_quads(
	SEXP a,
	SEXP b,
	SEXP c,
	SEXP d,
	const int **a_p,
	const int **b_p,
	const int **c_p,
	const int **d_p,
	const char *a_argname,
	const char *b_argname,
	const char *c_argname,
	const char *d_argname
);

SEXP Integer_sorted4(
	SEXP a,
	SEXP b,
	SEXP c,
	SEXP d,
	SEXP decreasing,
	SEXP strictly
);

SEXP Integer_order4(
	SEXP a,
	SEXP b,
	SEXP c,
	SEXP d,
	SEXP decreasing,
	SEXP use_radix
);

SEXP Integer_match4_quick(
	SEXP a1,
	SEXP b1,
	SEXP c1,
	SEXP d1,
	SEXP a2,
	SEXP b2,
	SEXP c2,
	SEXP d2,
	SEXP nomatch
);

SEXP Integer_selfmatch4_quick(
	SEXP a,
	SEXP b,
	SEXP c,
	SEXP d
);

SEXP Integer_match4_hash(
	SEXP a1,
	SEXP b1,
	SEXP c1,
	SEXP d1,
	SEXP a2,
	SEXP b2,
	SEXP c2,
	SEXP d2,
	SEXP nomatch
);

SEXP Integer_selfmatch4_hash(
	SEXP a,
	SEXP b,
	SEXP c,
	SEXP d
);

SEXP Integer_tabulate2(
	SEXP x,
	SEXP nbins,
	SEXP weight,
	SEXP strict
);

SEXP Integer_explode_bits(
	SEXP x,
	SEXP bitpos
);

SEXP Integer_sorted_merge(
	SEXP x,
	SEXP y
);

SEXP _find_interv_and_start_from_width(
	const int *x,
	int x_len,
	const int *width,
	int width_len
);

SEXP findIntervalAndStartFromWidth(
	SEXP x,
	SEXP vec
);


/* character_utils.c */

SEXP unstrsplit_list(SEXP x, SEXP sep);

SEXP safe_strexplode(SEXP s);

SEXP svn_time();


/* raw_utils.c */

SEXP _extract_bytes_by_positions(
	const char *x,
	int x_len,
	const int *pos,
	int npos,
	int collapse,
	SEXP lkup
);

SEXP _extract_bytes_by_ranges(
	const char *x,
	int x_len,
	const int *start,
	const int *width,
	int nranges,
	int collapse,
	SEXP lkup
);

SEXP C_extract_character_from_raw_by_positions(
	SEXP x,
	SEXP pos,
	SEXP collapse,
	SEXP lkup
);

SEXP C_extract_character_from_raw_by_ranges(
	SEXP x,
	SEXP start,
	SEXP width,
	SEXP collapse,
	SEXP lkup
);


/* eval_utils.c */

SEXP top_prenv(SEXP nm, SEXP env);

SEXP top_prenv_dots(SEXP env);


/* map_ranges_to_runs.c */

const char *_simple_range_mapper(
	const int *run_lengths,
	int nrun,
	int range_start,
	int range_end,
	int *mapped_range_offset,
	int *mapped_range_span,
	int *mapped_range_Ltrim,
	int *mapped_range_Rtrim
);

const char *_simple_position_mapper(
	const int *run_lengths,
	int nrun,
	int pos,
	int *mapped_pos
);

const char *_ranges_mapper(
	const int *run_lengths,
	int nrun,
	const int *start,
	const int *width,
	int nranges,
	int *mapped_range_offset,
	int *mapped_range_span,
	int *mapped_range_Ltrim,
	int *mapped_range_Rtrim,
	int method
);

const char *_positions_mapper(
	const int *run_lengths,
	int nrun,
	const int *pos,
	int npos,
	int *mapped_pos,
	int method
);

SEXP map_ranges(
	SEXP run_lengths,
	SEXP start,
	SEXP width,
	SEXP method
);

SEXP map_positions(
	SEXP run_lengths,
	SEXP pos,
	SEXP method
);


/* Hits_class.c */

SEXP _new_Hits(
	const char *Class,
	int *from,
	const int *to,
	int nhit,
	int nLnode,
	int nRnode,
	int already_sorted
);

SEXP Hits_new(
	SEXP Class,
	SEXP from,
	SEXP to,
	SEXP nLnode,
	SEXP nRnode,
	SEXP revmap_envir
);

int _get_select_mode(SEXP select);

SEXP select_hits(
	SEXP from,
	SEXP to,
	SEXP nLnode,
	SEXP nRnode,
	SEXP select,
	SEXP nodup
);

SEXP make_all_group_inner_hits(
	SEXP group_sizes,
	SEXP hit_type
);


/* Rle_class.c */

SEXP Rle_length(SEXP x);

SEXP Rle_valid(SEXP x);

SEXP _construct_logical_Rle(
	R_xlen_t nrun_in,
	const int *values_in,
	const void *lengths_in,
	int lengths_in_is_L
);

SEXP _construct_integer_Rle(
	R_xlen_t nrun_in,
	const int *values_in,
	const void *lengths_in,
	int lengths_in_is_L
);

SEXP _construct_numeric_Rle(
	R_xlen_t nrun_in,
	const double *values_in,
	const void *lengths_in,
	int lengths_in_is_L
);

SEXP _construct_complex_Rle(
	R_xlen_t nrun_in,
	const Rcomplex *values_in,
	const void *lengths_in,
	int lengths_in_is_L
);

SEXP _construct_character_Rle(
	SEXP values_in,
	const void *lengths_in,
	int lengths_in_is_L
);

SEXP _construct_raw_Rle(
	R_xlen_t nrun_in,
	const Rbyte *values_in,
	const void *lengths_in,
	int lengths_in_is_L
);

SEXP _construct_Rle(
	SEXP values_in,
	const void *lengths_in,
	int lengths_in_is_L
);

SEXP Rle_constructor(
	SEXP values_in,
	SEXP lengths_in
);

SEXP Rle_start(SEXP x);

SEXP Rle_end(SEXP x);

SEXP _subset_Rle_by_ranges(
	SEXP x,
	const int *start,
	const int *width,
	int nranges,
	int method,
	int as_list
);

SEXP _subset_Rle_by_positions(
	SEXP x,
	const int *pos,
	int npos,
	int method
);

SEXP Rle_extract_range(
	SEXP x,
	SEXP start,
	SEXP end
);

SEXP Rle_extract_ranges(
	SEXP x,
	SEXP start,
	SEXP width,
	SEXP method,
	SEXP as_list
);

SEXP Rle_extract_positions(
	SEXP x,
	SEXP pos,
	SEXP method
);

SEXP Rle_getStartEndRunAndOffset(
	SEXP x,
	SEXP start,
	SEXP end
);

SEXP Rle_window_aslist(
	SEXP x,
	SEXP runStart,
	SEXP runEnd,
	SEXP offsetStart,
	SEXP offsetEnd
);


/* Rle_utils.c */

SEXP Rle_runsum(
	SEXP x,
	SEXP k,
	SEXP na_rm
);

SEXP Rle_runwtsum(
	SEXP x,
	SEXP k,
	SEXP wt,
	SEXP na_rm
);

SEXP Rle_runq(
	SEXP x,
	SEXP k,
	SEXP which,
	SEXP na_rm
);


/* List_class.c */

const char *_get_List_elementType(SEXP x);

void _set_List_elementType(
	SEXP x,
	const char *type
);


/* SimpleList_class.c */

SEXP _new_SimpleList(
	const char *classname,
	SEXP listData
);


/* DataFrame_class.c */

SEXP _new_DataFrame(
	const char *classname,
	SEXP vars,
	SEXP rownames,
	SEXP nrows
);

