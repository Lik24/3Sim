/*
 * vibor_t_terminate.c
 *
 * Code generation for function 'vibor_t_terminate'
 *
 * C source code generated on: Tue Aug 12 14:07:28 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "vibor_t.h"
#include "vibor_t_terminate.h"

/* Function Definitions */
void vibor_t_atexit(void)
{
  emlrtCreateRootTLS(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1);
  emlrtEnterRtStackR2012b(emlrtRootTLSGlobal);
  emlrtLeaveRtStackR2012b(emlrtRootTLSGlobal);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

void vibor_t_terminate(void)
{
  emlrtLeaveRtStackR2012b(emlrtRootTLSGlobal);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

/* End of code generation (vibor_t_terminate.c) */
