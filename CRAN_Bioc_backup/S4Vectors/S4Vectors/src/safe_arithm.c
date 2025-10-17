/****************************************************************************
 * Safe signed integer arithmetic                                           *
 * ------------------------------                                           *
 * TODO: Extend to support safe double arithmetic when the need arises.     *
 ****************************************************************************/
#include "S4Vectors.h"

#include <limits.h>  /* for INT_MAX, INT_MIN, LLONG_MAX, and LLONG_MIN */
#include <ctype.h>   /* for isdigit() and isspace() */

static int ovflow_flag;

void _reset_ovflow_flag()
{
	ovflow_flag = 0;
	return;
}

int _get_ovflow_flag()
{
	return ovflow_flag;
}


/****************************************************************************
 * Safe arithmetic on int values
 *
 * Reference:
 *   The CERT C Secure Coding Standard
 *     Rule INT32-C. Ensure that operations on signed integers do not result
 *     in overflow
 */

int _safe_int_add(int x, int y)
{
	if (x == NA_INTEGER || y == NA_INTEGER)
		return NA_INTEGER;
	if ((y > 0 && x > INT_MAX - y) ||
	    (y < 0 && x < INT_MIN - y))
	{
		ovflow_flag = 1;
		return NA_INTEGER;
	}
	return x + y;
}

int _safe_int_subtract(int x, int y)
{
	if (x == NA_INTEGER || y == NA_INTEGER)
		return NA_INTEGER;
	if ((y < 0 && x > INT_MAX + y) ||
	    (y > 0 && x < INT_MIN + y))
	{
		ovflow_flag = 1;
		return NA_INTEGER;
	}
	return x - y;
}

int _safe_int_mult(int x, int y)
{
	if (x == NA_INTEGER || y == NA_INTEGER)
		return NA_INTEGER;
	if (x > 0) { /* x is positive */
		if (y > 0) { /* x and y are positive */
			if (x > (INT_MAX / y)) {
				ovflow_flag = 1;
				return NA_INTEGER;
			}
		} else { /* x is positive, y is non-positive */
			if (y < (INT_MIN / x)) {
				ovflow_flag = 1;
				return NA_INTEGER;
			}
		}
	} else { /* x is non-positive */
		if (y > 0) { /* x is non-positive, y is positive */
			if (x < (INT_MIN / y)) {
				ovflow_flag = 1;
				return NA_INTEGER;
			}
	  	} else { /* x and y are non-positive */
			if ((x != 0) && (y < (INT_MAX / x))) {
				ovflow_flag = 1;
				return NA_INTEGER;
			}
		}
	}
	return x * y;
}


/****************************************************************************
 * _as_int()
 *
 * Turn string pointed by 'val' into an int. The string has no terminating
 * null byte ('\0') and must have the following format:
 *     ^[[:space:]]*[+-]?[[:digit:]]+[[:space:]]*$
 * Return NA_INTEGER if the string is malformed or if it represents an integer
 * value that cannot be represented by an int (int overflow).
 * TODO: Maybe implement this on top of strtol(). Would be much simpler but
 * would it be equivalent? Also would it be as fast? See how as_double() in
 * rtracklayer/src/readGFF.c is implemented on top of strtod().
 */
#define LEADING_SPACE 0
#define NUMBER 1
#define TRAILING_SPACE 2
int _as_int(const char *val, int val_len)
{
	int n, ndigit, sign, state, i;
	char c;

	n = ndigit = 0;
	sign = 1;
	state = LEADING_SPACE;
	for (i = 0; i < val_len; i++) {
		c = val[i];
		if (isdigit(c)) {
			if (state == TRAILING_SPACE)
				return NA_INTEGER;  /* malformed string */
			state = NUMBER;
			ndigit++;
			n = _safe_int_mult(n, 10);
			n = _safe_int_add(n, c - '0');
			if (n == NA_INTEGER)
				return NA_INTEGER;  /* int overflow */
			continue;
		}
		if (c == '+' || c == '-') {
			if (state != LEADING_SPACE)
				return NA_INTEGER;  /* malformed string */
			state = NUMBER;
			if (c == '-')
				sign = -1;
			continue;
		}
		if (!isspace(c))
			return NA_INTEGER;  /* malformed string */
		if (state == NUMBER) {
			if (ndigit == 0)
				return NA_INTEGER;  /* malformed string */
			state = TRAILING_SPACE;
		}
	}
	if (ndigit == 0)
		return NA_INTEGER;  /* malformed string */
	if (sign == -1)
		n = -n;
	return n;
}


/****************************************************************************
 * Safe arithmetic on long long int values
 */

long long int _safe_llint_add(long long int x, long long int y)
{
	if (x == NA_LLINT || y == NA_LLINT)
		return NA_LLINT;
	if ((y > 0LL && x > LLONG_MAX - y) ||
	    (y < 0LL && x < LLONG_MIN - y))
	{
		ovflow_flag = 1;
		return NA_LLINT;
	}
	return x + y;
}

long long int _safe_llint_subtract(long long int x, long long int y)
{
	if (x == NA_LLINT || y == NA_LLINT)
		return NA_LLINT;
	if ((y < 0LL && x > LLONG_MAX + y) ||
	    (y > 0LL && x < LLONG_MIN + y))
	{
		ovflow_flag = 1;
		return NA_LLINT;
	}
	return x - y;
}

long long int _safe_llint_mult(long long int x, long long int y)
{
	if (x == NA_LLINT || y == NA_LLINT)
		return NA_LLINT;
	if (x > 0LL) { /* x is positive */
		if (y > 0LL) { /* x and y are positive */
			if (x > (LLONG_MAX / y)) {
				ovflow_flag = 1;
				return NA_LLINT;
			}
		} else { /* x is positive, y is non-positive */
			if (y < (LLONG_MIN / x)) {
				ovflow_flag = 1;
				return NA_LLINT;
			}
		}
	} else { /* x is non-positive */
		if (y > 0LL) { /* x is non-positive, y is positive */
			if (x < (LLONG_MIN / y)) {
				ovflow_flag = 1;
				return NA_LLINT;
			}
	  	} else { /* x and y are non-positive */
			if ((x != 0LL) && (y < (LLONG_MAX / x))) {
				ovflow_flag = 1;
				return NA_LLINT;
			}
		}
	}
	return x * y;
}

