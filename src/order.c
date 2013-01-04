/**
 * Routines that return a permutation of indeces
 * to put an array in sorted order. These are similar
 * to R's order function.
 *
 * David Faden, dfaden@gmail.com
 * August 21, 2005
 */

#include <assert.h>
#include "order.h"
#include <stdlib.h>
#include <string.h>

typedef struct {
  size_t index;

  /**
   * Pointer to datum.
   */
  const void* datum;

  /**
   * Pointer to comparison function.
   */
  ComparisonFunc compare;
} Pair;

int comparePairs(const void* v1, const void* v2)
{
  const Pair* p1 = v1;
  const Pair* p2 = v2;
  assert(p1->compare == p2->compare);
  return p1->compare(p1->datum, p2->datum);
}

/**
 * The caller is responsible for freeing the returned object.
 */
size_t* order(const void* base, size_t numElements, size_t size,
	      ComparisonFunc compare)
{
  Pair* pairs = malloc(sizeof(Pair) * numElements);
  size_t* indeces = malloc(sizeof(size_t) * numElements);
  size_t i;
  const char* p = base;

  for (i = 0; i < numElements; ++i, p += size) {
    pairs[i].index = i;
    pairs[i].datum = p;
    pairs[i].compare = compare;
  }

  qsort(pairs, numElements, sizeof(Pair), comparePairs);

  for (i = 0; i < numElements; ++i) {
    indeces[i] = pairs[i].index;
  }

  free(pairs);

  return indeces;
}

#define RETURN_CMP(a,b) if (a < b) { return -1; } \
  else if (a == b) { return 0; } \
  else { return 1; }

int intCompare(const void* v1, const void* v2)
{
  int i1 = *((const int*)v1);
  int i2 = *((const int*)v2);
  RETURN_CMP(i1, i2);
}

int doubleCompare(const void* v1, const void* v2)
{
  double d1 = *((const double*)v1);
  double d2 = *((const double*)v2);
  RETURN_CMP(d1, d2);
}

int stringCompare(const void* v1, const void* v2)
{
  return strcmp(*(const char**)v1, *(const char**)v2);
}

size_t* orderInt(const int* base, size_t numElements)
{
  return order(base, numElements, sizeof(int), intCompare);
}

size_t* orderDouble(const double* base, size_t numElements)
{
  return order(base, numElements, sizeof(double), doubleCompare); 
}

size_t* orderString(char* const* base, size_t numElements)
{
  return order(base, numElements, sizeof(const char*), stringCompare);
}


