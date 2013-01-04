/**
 * Routines that return a permutation of indeces
 * to put an array in sorted order. These are similar
 * to R's order function.
 *
 * David Faden, dfaden@gmail.com
 * August 21, 2005
 */

#ifndef __ORDER_H__
#define __ORDER_H__

#include <stddef.h>

size_t* orderInt(const int* base, size_t numElements);

size_t* orderDouble(const double* base, size_t numElements);

size_t* orderString(char* const* base, size_t numElements);

/**
 * Returns a value < 0 if v1 is "less than" v2.
 * Returns 0 if v1 "equals" v2.
 * Returns a value > 0 if v1 is "greater than" v2.
 */
typedef int (*ComparisonFunc)(const void* v1, const void* v2);

/**
 * Returns a permutation of indeces for elements
 * of base that would bring it into ascending sorted
 * order.
 *
 * Notice that the parameters have the same meanings
 * as for stdlib's qsort.
 *
 * The caller is responsible for freeing the returned object.
 *
 * base -- pointer to the first element in the array to be ordered.
 * numElements -- number of elements in the array.
 * size -- the size in bytes of element in array, as may be determined
 * with sizeof.
 * compare -- comparison function.
 */
size_t* order(const void* base, size_t numElements, size_t size,
	      ComparisonFunc compare);


/**
 * This macro declares a function named order_TYPE
 * taking a const TYPE* as its first argument
 * and a size_t as its second -- the array
 * and the length of the array to be sorted
 * respectively.
 */
#define MAKE_ORDER_FUNC(TYPE) size_t* order_ ## TYPE (const TYPE* a, size_t len) \
  { \
    int compare_ ## TYPE (const void* v1, const void* v2) \
    { \
      TYPE t1 = *((const TYPE*)v1);		\
      TYPE t2 = *((const TYPE*)v2);		\
      if (t1 < t2) { return -1; }		\
      else if (t1 == t2) { return 0; }		\
      else { return 1; }				\
    }								\
    return order(a, len, sizeof( TYPE ), compare_ ## TYPE );	\
}


#endif
