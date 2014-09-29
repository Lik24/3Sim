/*
 * rdivide.c
 *
 * Code generation for function 'rdivide'
 *
 * C source code generated on: Mon Jul 21 16:57:34 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "vibor_t.h"
#include "rdivide.h"
#include "vibor_t_emxutil.h"
#include "vibor_t_mexutil.h"

/* Variable Definitions */
static emlrtRSInfo s_emlrtRSI = { 13, "rdivide",
  "C:/Program Files/MATLAB/R2013a/toolbox/eml/lib/matlab/ops/rdivide.m" };

static emlrtMCInfo h_emlrtMCI = { 14, 5, "rdivide",
  "C:/Program Files/MATLAB/R2013a/toolbox/eml/lib/matlab/ops/rdivide.m" };

static emlrtMCInfo i_emlrtMCI = { 13, 15, "rdivide",
  "C:/Program Files/MATLAB/R2013a/toolbox/eml/lib/matlab/ops/rdivide.m" };

static emlrtRTEInfo c_emlrtRTEI = { 1, 14, "rdivide",
  "C:/Program Files/MATLAB/R2013a/toolbox/eml/lib/matlab/ops/rdivide.m" };

/* Function Definitions */
void rdivide(const emxArray_real_T *x, const emxArray_real_T *y, emxArray_real_T
             *z)
{
  boolean_T p;
  boolean_T b_p;
  int32_T i;
  boolean_T exitg1;
  int32_T b_x[2];
  int32_T b_y[2];
  const mxArray *c_y;
  static const int32_T iv2[2] = { 1, 15 };

  const mxArray *m2;
  char_T cv4[15];
  static const char_T cv5[15] = { 'M', 'A', 'T', 'L', 'A', 'B', ':', 'd', 'i',
    'm', 'a', 'g', 'r', 'e', 'e' };

  int32_T loop_ub;
  p = FALSE;
  b_p = TRUE;
  i = 0;
  exitg1 = FALSE;
  while ((exitg1 == FALSE) && (i < 2)) {
    b_x[0] = x->size[0];
    b_x[1] = 1;
    b_y[0] = y->size[0];
    b_y[1] = 1;
    if (!(b_x[i] == b_y[i])) {
      b_p = FALSE;
      exitg1 = TRUE;
    } else {
      i++;
    }
  }

  if (!b_p) {
  } else {
    p = TRUE;
  }

  if (p) {
  } else {
    emlrtPushRtStackR2012b(&s_emlrtRSI, emlrtRootTLSGlobal);
    c_y = NULL;
    m2 = mxCreateCharArray(2, iv2);
    for (i = 0; i < 15; i++) {
      cv4[i] = cv5[i];
    }

    emlrtInitCharArrayR2013a(emlrtRootTLSGlobal, 15, m2, cv4);
    emlrtAssign(&c_y, m2);
    error(message(c_y, &h_emlrtMCI), &i_emlrtMCI);
    emlrtPopRtStackR2012b(&s_emlrtRSI, emlrtRootTLSGlobal);
  }

  i = z->size[0];
  z->size[0] = x->size[0];
  emxEnsureCapacity((emxArray__common *)z, i, (int32_T)sizeof(real_T),
                    &c_emlrtRTEI);
  loop_ub = x->size[0];
  for (i = 0; i < loop_ub; i++) {
    z->data[i] = x->data[i] / y->data[i];
  }
}

/* End of code generation (rdivide.c) */
