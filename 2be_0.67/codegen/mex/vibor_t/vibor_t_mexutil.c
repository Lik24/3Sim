/*
 * vibor_t_mexutil.c
 *
 * Code generation for function 'vibor_t_mexutil'
 *
 * C source code generated on: Tue Aug 12 14:07:28 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "vibor_t.h"
#include "vibor_t_mexutil.h"

/* Function Definitions */
void error(const mxArray *b, emlrtMCInfo *location)
{
  const mxArray *pArray;
  pArray = b;
  emlrtCallMATLABR2012b(emlrtRootTLSGlobal, 0, NULL, 1, &pArray, "error", TRUE,
                        location);
}

const mxArray *message(const mxArray *b, emlrtMCInfo *location)
{
  const mxArray *pArray;
  const mxArray *m6;
  pArray = b;
  return emlrtCallMATLABR2012b(emlrtRootTLSGlobal, 1, &m6, 1, &pArray, "message",
    TRUE, location);
}

/* End of code generation (vibor_t_mexutil.c) */
