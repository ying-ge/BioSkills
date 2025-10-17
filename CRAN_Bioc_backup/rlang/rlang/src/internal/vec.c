#include <rlang.h>
#include "vec.h"
#include "decl/vec-decl.h"


bool r_is_atomic(r_obj* x, r_ssize n) {
  switch(r_typeof(x)) {
  case R_TYPE_logical:
  case R_TYPE_integer:
  case R_TYPE_double:
  case R_TYPE_complex:
  case R_TYPE_character:
  case RAWSXP:
    return _r_has_correct_length(x, n);
  default:
    return false;
  }
}

bool r_is_vector(r_obj* x, r_ssize n) {
  switch(r_typeof(x)) {
  case R_TYPE_logical:
  case R_TYPE_integer:
  case R_TYPE_double:
  case R_TYPE_complex:
  case R_TYPE_character:
  case RAWSXP:
  case VECSXP:
    return _r_has_correct_length(x, n);
  default:
    return false;
  }
}

bool r_is_logical(r_obj* x, r_ssize n) {
  return r_typeof(x) == R_TYPE_logical && _r_has_correct_length(x, n);
}

bool r_is_integer(r_obj* x, r_ssize n, int finite) {
  if (r_typeof(x) != R_TYPE_integer || !_r_has_correct_length(x, n)) {
    return false;
  }
  if (finite >= 0 && (bool) finite != _r_is_finite(x)) {
    return false;
  }
  return true;
}
bool r_is_double(r_obj* x, r_ssize n, int finite) {
  return _r_is_double(x, n, finite);
}
bool r_is_complex(r_obj* x, r_ssize n, int finite) {
  return _r_is_complex(x, n, finite);
}

bool r_is_integerish(r_obj* x, r_ssize n, int finite) {
  if (r_typeof(x) == R_TYPE_integer) {
    return r_is_integer(x, n, finite);
  }
  if (r_typeof(x) != R_TYPE_double || !_r_has_correct_length(x, n)) {
    return false;
  }

  r_ssize actual_n = r_length(x);
  const double* p_x = r_dbl_cbegin(x);
  bool actual_finite = true;

  for (r_ssize i = 0; i < actual_n; ++i) {
    double elt = p_x[i];

    if (!isfinite(elt)) {
      actual_finite = false;
      continue;
    }

    if (!r_dbl_is_whole(elt)) {
      return false;
    }
  }

  if (finite >= 0 && actual_finite != (bool) finite) {
    return false;
  }

  return true;
}

bool is_character(r_obj* x,
                  r_ssize n,
                  enum option_bool missing,
                  enum option_bool empty) {
  if (r_typeof(x) != R_TYPE_character) {
    return false;
  }
  if (!_r_has_correct_length(x, n)) {
    return false;
  }

  bool has_missing = missing != OPTION_BOOL_null;
  bool has_empty = empty != OPTION_BOOL_null;

  if (!has_missing && !has_empty) {
    return true;
  }
  if (missing == OPTION_BOOL_true && empty == OPTION_BOOL_true) {
    r_abort("Exactly one of `missing` and `empty` can be `TRUE`.");
  }

  n = r_length(x);
  r_obj* const * v_x = r_chr_cbegin(x);

  // Could we inspect ALTREP properties for the `missing` case?
  if (!list_match(v_x, n, r_strs.na, missing)) {
    return false;
  }
  if (!list_match(v_x, n, r_strs.empty, empty)) {
    return false;
  }

  return true;
}

static
bool list_match(r_obj* const * v_x,
                r_ssize n,
                r_obj* value,
                enum option_bool match) {
  switch (match) {
  case OPTION_BOOL_null:
    return true;
  case OPTION_BOOL_true:
    for (r_ssize i = 0; i < n; ++i) {
      if (v_x[i] != value) {
        return false;
      }
    }
    return true;
  case OPTION_BOOL_false:
    for (r_ssize i = 0; i < n; ++i) {
      if (v_x[i] == value) {
        return false;
      }
    }
    return true;
  default:
    r_stop_unreachable();
  }
}

bool r_is_raw(r_obj* x, r_ssize n) {
  return r_typeof(x) == R_TYPE_raw && _r_has_correct_length(x, n);
}

r_ssize validate_n(r_obj* n) {
  if (n == r_null) {
    return -1;
  }

  switch (r_typeof(n)) {
  case R_TYPE_integer:
  case R_TYPE_double:
    if (r_length(n) == 1) {
      break;
    }
    // fallthrough
  default:
    r_abort("`n` must be NULL or a scalar integer");
  }

  return r_arg_as_ssize(n, "n");
}


// Coercion ----------------------------------------------------------

static
r_obj* vec_coercer(r_obj* to) {
  switch (r_typeof(to)) {
  case R_TYPE_logical: return rlang_ns_get("legacy_as_logical");
  case R_TYPE_integer: return rlang_ns_get("legacy_as_integer");
  case R_TYPE_double: return rlang_ns_get("legacy_as_double");
  case R_TYPE_complex: return rlang_ns_get("legacy_as_complex");
  case R_TYPE_character: return rlang_ns_get("legacy_as_character");
  case RAWSXP: return rlang_ns_get("legacy_as_raw");
  default: r_abort("No coercion implemented for `%s`", Rf_type2str(r_typeof(to)));
  }
}

void r_vec_poke_coerce_n(r_obj* x, r_ssize offset,
                         r_obj* y, r_ssize from, r_ssize n) {
  if (r_typeof(y) == r_typeof(x)) {
    r_vec_poke_n(x, offset, y, from, n);
    return ;
  }
  if (r_is_object(y)) {
    r_abort("Can't splice S3 objects");
  }

  // FIXME: This callbacks to rlang R coercers with an extra copy.
  r_obj* coercer = vec_coercer(x);
  r_obj* call = KEEP(Rf_lang2(coercer, y));
  r_obj* coerced = KEEP(r_eval(call, R_BaseEnv));

  r_vec_poke_n(x, offset, coerced, from, n);
  FREE(2);
}

void r_vec_poke_coerce_range(r_obj* x, r_ssize offset,
                             r_obj* y, r_ssize from, r_ssize to) {
  r_vec_poke_coerce_n(x, offset, y, from, to - from + 1);
}
