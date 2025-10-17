#include "S4Vectors.h"

#define CALLMETHOD_DEF(fun, numArgs) {#fun, (DL_FUNC) &fun, numArgs}

#define REGISTER_CCALLABLE(fun) \
	R_RegisterCCallable("S4Vectors", #fun, (DL_FUNC) &fun)


static const R_CallMethodDef callMethods[] = {

/* sort_utils.c */
	CALLMETHOD_DEF(test_sort_ushort_array, 2),

/* AEbufs.c */
	CALLMETHOD_DEF(AEbufs_use_malloc, 1),
	CALLMETHOD_DEF(AEbufs_free, 0),

/* anyMissing.c */
	CALLMETHOD_DEF(anyMissing, 1),

/* LLint_class.c */
	CALLMETHOD_DEF(make_RAW_from_NA_LLINT, 0),
	CALLMETHOD_DEF(new_LLint_from_LOGICAL, 1),
	CALLMETHOD_DEF(new_LLint_from_INTEGER, 1),
	CALLMETHOD_DEF(new_LLint_from_NUMERIC, 1),
	CALLMETHOD_DEF(new_LLint_from_CHARACTER, 1),
	CALLMETHOD_DEF(new_LOGICAL_from_LLint, 1),
	CALLMETHOD_DEF(new_INTEGER_from_LLint, 1),
	CALLMETHOD_DEF(new_NUMERIC_from_LLint, 1),
	CALLMETHOD_DEF(new_CHARACTER_from_LLint, 1),

/* subsetting_utils.c */
	CALLMETHOD_DEF(vector_OR_factor_extract_positions, 2),
	CALLMETHOD_DEF(vector_OR_factor_extract_ranges, 3),

/* vector_utils.c */
	CALLMETHOD_DEF(sapply_NROW, 1),

/* logical_utils.c */
	CALLMETHOD_DEF(logical_sum, 2),
	CALLMETHOD_DEF(logical2_sum, 2),

/* integer_utils.c */
	CALLMETHOD_DEF(to_list_of_ints, 2),
	CALLMETHOD_DEF(Integer_any_missing_or_outside, 3),
	CALLMETHOD_DEF(Integer_diff_with_0, 1),
	CALLMETHOD_DEF(Integer_diff_with_last, 2),
	CALLMETHOD_DEF(Integer_order, 3),
	CALLMETHOD_DEF(Integer_pcompare2, 4),
	CALLMETHOD_DEF(Integer_sorted2, 4),
	CALLMETHOD_DEF(Integer_order2, 4),
	CALLMETHOD_DEF(Integer_match2_quick, 5),
	CALLMETHOD_DEF(Integer_selfmatch2_quick, 2),
	CALLMETHOD_DEF(Integer_match2_hash, 5),
	CALLMETHOD_DEF(Integer_selfmatch2_hash, 2),
	CALLMETHOD_DEF(Integer_sorted4, 6),
	CALLMETHOD_DEF(Integer_order4, 6),
	CALLMETHOD_DEF(Integer_match4_quick, 9),
	CALLMETHOD_DEF(Integer_selfmatch4_quick, 4),
	CALLMETHOD_DEF(Integer_match4_hash, 9),
	CALLMETHOD_DEF(Integer_selfmatch4_hash, 4),
	CALLMETHOD_DEF(Integer_tabulate2, 4),
	CALLMETHOD_DEF(Integer_explode_bits, 2),
	CALLMETHOD_DEF(Integer_sorted_merge, 2),
	CALLMETHOD_DEF(findIntervalAndStartFromWidth, 2),

/* character_utils.c */
	CALLMETHOD_DEF(unstrsplit_list, 2),
	CALLMETHOD_DEF(safe_strexplode, 1),
	CALLMETHOD_DEF(svn_time, 0),

/* raw_utils.c */
	CALLMETHOD_DEF(C_extract_character_from_raw_by_positions, 4),
	CALLMETHOD_DEF(C_extract_character_from_raw_by_ranges, 5),

/* eval_utils.c */
	CALLMETHOD_DEF(top_prenv, 2),
	CALLMETHOD_DEF(top_prenv_dots, 1),

/* map_ranges_to_runs.c */
	CALLMETHOD_DEF(map_ranges, 4),
	CALLMETHOD_DEF(map_positions, 3),

/* Hits_class.c */
	CALLMETHOD_DEF(Hits_new, 6),
	CALLMETHOD_DEF(select_hits, 6),
	CALLMETHOD_DEF(make_all_group_inner_hits, 2),

/* Rle_class.c */
	CALLMETHOD_DEF(Rle_length, 1),
	CALLMETHOD_DEF(Rle_valid, 1),
	CALLMETHOD_DEF(Rle_constructor, 2),
	CALLMETHOD_DEF(Rle_start, 1),
	CALLMETHOD_DEF(Rle_end, 1),
	CALLMETHOD_DEF(Rle_extract_range, 3),
	CALLMETHOD_DEF(Rle_extract_ranges, 5),
	CALLMETHOD_DEF(Rle_extract_positions, 3),
	CALLMETHOD_DEF(Rle_getStartEndRunAndOffset, 3),
	CALLMETHOD_DEF(Rle_window_aslist, 5),

/* Rle_utils.c */
	CALLMETHOD_DEF(Rle_runsum, 3),
	CALLMETHOD_DEF(Rle_runwtsum, 4),
	CALLMETHOD_DEF(Rle_runq, 4),

	{NULL, NULL, 0}
};


void R_init_S4Vectors(DllInfo *info)
{
	R_registerRoutines(info, NULL, callMethods, NULL, NULL);

/* safe_arithm.c */
	REGISTER_CCALLABLE(_reset_ovflow_flag);
	REGISTER_CCALLABLE(_get_ovflow_flag);
	REGISTER_CCALLABLE(_safe_int_add);
	REGISTER_CCALLABLE(_safe_int_mult);
	REGISTER_CCALLABLE(_as_int);
	REGISTER_CCALLABLE(_safe_llint_add);
	REGISTER_CCALLABLE(_safe_llint_mult);

/* sort_utils.c */
	REGISTER_CCALLABLE(_sort_ints);
	REGISTER_CCALLABLE(_get_order_of_int_array);
	REGISTER_CCALLABLE(_sort_int_array);
	REGISTER_CCALLABLE(_get_order_of_int_pairs);
	REGISTER_CCALLABLE(_sort_int_pairs);
	REGISTER_CCALLABLE(_get_matches_of_ordered_int_pairs);
	REGISTER_CCALLABLE(_get_order_of_int_quads);
	REGISTER_CCALLABLE(_get_matches_of_ordered_int_quads);

/* hash_utils.c */
	REGISTER_CCALLABLE(_new_htab);
	REGISTER_CCALLABLE(_get_hbucket_val);
	REGISTER_CCALLABLE(_set_hbucket_val);

/* AEbufs.c */
	REGISTER_CCALLABLE(_increase_buflength);
	REGISTER_CCALLABLE(_IntAE_get_nelt);
	REGISTER_CCALLABLE(_IntAE_set_nelt);
	REGISTER_CCALLABLE(_IntAE_set_val);
	REGISTER_CCALLABLE(_IntAE_extend);
	REGISTER_CCALLABLE(_IntAE_insert_at);
	REGISTER_CCALLABLE(_new_IntAE);
	REGISTER_CCALLABLE(_IntAE_append);
	REGISTER_CCALLABLE(_IntAE_delete_at);
	REGISTER_CCALLABLE(_IntAE_shift);
	REGISTER_CCALLABLE(_IntAE_sum_and_shift);
	REGISTER_CCALLABLE(_IntAE_qsort);
	REGISTER_CCALLABLE(_IntAE_uniq);
	REGISTER_CCALLABLE(_new_INTEGER_from_IntAE);
	REGISTER_CCALLABLE(_new_LOGICAL_from_IntAE);
	REGISTER_CCALLABLE(_new_IntAE_from_INTEGER);
	REGISTER_CCALLABLE(_new_IntAE_from_CHARACTER);
	REGISTER_CCALLABLE(_IntAEAE_get_nelt);
	REGISTER_CCALLABLE(_IntAEAE_set_nelt);
	REGISTER_CCALLABLE(_IntAEAE_extend);
	REGISTER_CCALLABLE(_IntAEAE_insert_at);
	REGISTER_CCALLABLE(_new_IntAEAE);
	REGISTER_CCALLABLE(_IntAEAE_pappend);
	REGISTER_CCALLABLE(_IntAEAE_shift);
	REGISTER_CCALLABLE(_IntAEAE_sum_and_shift);
	REGISTER_CCALLABLE(_new_LIST_from_IntAEAE);
	REGISTER_CCALLABLE(_new_IntAEAE_from_LIST);
	REGISTER_CCALLABLE(_IntAEAE_toEnvir);
	REGISTER_CCALLABLE(_IntPairAE_get_nelt);
	REGISTER_CCALLABLE(_IntPairAE_set_nelt);
	REGISTER_CCALLABLE(_IntPairAE_extend);
	REGISTER_CCALLABLE(_IntPairAE_insert_at);
	REGISTER_CCALLABLE(_new_IntPairAE);
	REGISTER_CCALLABLE(_IntPairAEAE_get_nelt);
	REGISTER_CCALLABLE(_IntPairAEAE_set_nelt);
	REGISTER_CCALLABLE(_IntPairAEAE_extend);
	REGISTER_CCALLABLE(_IntPairAEAE_insert_at);
	REGISTER_CCALLABLE(_new_IntPairAEAE);
	REGISTER_CCALLABLE(_LLongAE_get_nelt);
	REGISTER_CCALLABLE(_LLongAE_set_nelt);
	REGISTER_CCALLABLE(_LLongAE_set_val);
	REGISTER_CCALLABLE(_LLongAE_extend);
	REGISTER_CCALLABLE(_LLongAE_insert_at);
	REGISTER_CCALLABLE(_new_LLongAE);
	REGISTER_CCALLABLE(_LLongAEAE_get_nelt);
	REGISTER_CCALLABLE(_LLongAEAE_set_nelt);
	REGISTER_CCALLABLE(_LLongAEAE_extend);
	REGISTER_CCALLABLE(_LLongAEAE_insert_at);
	REGISTER_CCALLABLE(_new_LLongAEAE);
	REGISTER_CCALLABLE(_DoubleAE_get_nelt);
	REGISTER_CCALLABLE(_DoubleAE_set_nelt);
	REGISTER_CCALLABLE(_DoubleAE_set_val);
	REGISTER_CCALLABLE(_DoubleAE_extend);
	REGISTER_CCALLABLE(_DoubleAE_insert_at);
	REGISTER_CCALLABLE(_new_DoubleAE);
	REGISTER_CCALLABLE(_DoubleAE_append);
	REGISTER_CCALLABLE(_DoubleAE_delete_at);
	REGISTER_CCALLABLE(_new_NUMERIC_from_DoubleAE);
	REGISTER_CCALLABLE(_new_DoubleAE_from_NUMERIC);
	REGISTER_CCALLABLE(_CharAE_get_nelt);
	REGISTER_CCALLABLE(_CharAE_set_nelt);
	REGISTER_CCALLABLE(_CharAE_extend);
	REGISTER_CCALLABLE(_CharAE_insert_at);
	REGISTER_CCALLABLE(_new_CharAE);
	REGISTER_CCALLABLE(_new_CharAE_from_string);
	REGISTER_CCALLABLE(_CharAE_append_string);
	REGISTER_CCALLABLE(_CharAE_delete_at);
	REGISTER_CCALLABLE(_new_CHARSXP_from_CharAE);
	REGISTER_CCALLABLE(_new_RAW_from_CharAE);
	REGISTER_CCALLABLE(_new_LOGICAL_from_CharAE);
	REGISTER_CCALLABLE(_CharAEAE_get_nelt);
	REGISTER_CCALLABLE(_CharAEAE_set_nelt);
	REGISTER_CCALLABLE(_CharAEAE_extend);
	REGISTER_CCALLABLE(_CharAEAE_insert_at);
	REGISTER_CCALLABLE(_new_CharAEAE);
	REGISTER_CCALLABLE(_CharAEAE_append_string);
	REGISTER_CCALLABLE(_new_CHARACTER_from_CharAEAE);

/* SEXP_utils.c */
	REGISTER_CCALLABLE(_get_classname);

/* LLint_class.c */
	REGISTER_CCALLABLE(_is_LLint);
	REGISTER_CCALLABLE(_get_LLint_length);
	REGISTER_CCALLABLE(_get_LLint_dataptr);
	REGISTER_CCALLABLE(_alloc_LLint);

/* subsetting_utils.c */
	REGISTER_CCALLABLE(_copy_vector_block);
	REGISTER_CCALLABLE(_copy_vector_positions);
	REGISTER_CCALLABLE(_copy_vector_ranges);

/* vector_utils.c */
	REGISTER_CCALLABLE(_vector_memcmp);
	REGISTER_CCALLABLE(_list_as_data_frame);

/* integer_utils.c */
	REGISTER_CCALLABLE(_check_integer_pairs);
	REGISTER_CCALLABLE(_find_interv_and_start_from_width);

/* raw_utils.c */
	REGISTER_CCALLABLE(_extract_bytes_by_positions);
	REGISTER_CCALLABLE(_extract_bytes_by_ranges);

/* Hits_class.c */
	REGISTER_CCALLABLE(_new_Hits);
	REGISTER_CCALLABLE(_get_select_mode);

/* Rle_class.c */
	REGISTER_CCALLABLE(_construct_logical_Rle);
	REGISTER_CCALLABLE(_construct_integer_Rle);
	REGISTER_CCALLABLE(_construct_numeric_Rle);
	REGISTER_CCALLABLE(_construct_complex_Rle);
	REGISTER_CCALLABLE(_construct_character_Rle);
	REGISTER_CCALLABLE(_construct_raw_Rle);
	REGISTER_CCALLABLE(_construct_Rle);

/* List_class.c */
	REGISTER_CCALLABLE(_get_List_elementType);
	REGISTER_CCALLABLE(_set_List_elementType);

/* SimpleList_class.c */
	REGISTER_CCALLABLE(_new_SimpleList);

/* DataFrame_class.c */
	REGISTER_CCALLABLE(_new_DataFrame);

	return;
}

