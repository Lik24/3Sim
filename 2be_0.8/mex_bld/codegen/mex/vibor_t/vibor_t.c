/*
 * vibor_t.c
 *
 * Code generation for function 'vibor_t'
 *
 * C source code generated on: Mon Jul 21 16:57:34 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "vibor_t.h"
#include "vibor_t_emxutil.h"
#include "abs.h"
#include "rdivide.h"
#include "vibor_t_mexutil.h"

/* Variable Definitions */
static emlrtRSInfo emlrtRSI = { 15, "vibor_t",
  "C:/Users/LIK/Documents/MATLAB/Tube_dev/2be_0.660/Sim_Lib/vibor_t.m" };

static emlrtRSInfo b_emlrtRSI = { 18, "vibor_t",
  "C:/Users/LIK/Documents/MATLAB/Tube_dev/2be_0.660/Sim_Lib/vibor_t.m" };

static emlrtRSInfo c_emlrtRSI = { 25, "vibor_t",
  "C:/Users/LIK/Documents/MATLAB/Tube_dev/2be_0.660/Sim_Lib/vibor_t.m" };

static emlrtRSInfo d_emlrtRSI = { 26, "vibor_t",
  "C:/Users/LIK/Documents/MATLAB/Tube_dev/2be_0.660/Sim_Lib/vibor_t.m" };

static emlrtRSInfo e_emlrtRSI = { 16, "min",
  "C:/Program Files/MATLAB/R2013a/toolbox/eml/lib/matlab/datafun/min.m" };

static emlrtRSInfo f_emlrtRSI = { 18, "eml_min_or_max",
  "C:/Program Files/MATLAB/R2013a/toolbox/eml/lib/matlab/eml/eml_min_or_max.m" };

static emlrtRSInfo g_emlrtRSI = { 38, "eml_min_or_max",
  "C:/Program Files/MATLAB/R2013a/toolbox/eml/lib/matlab/eml/eml_min_or_max.m" };

static emlrtRSInfo h_emlrtRSI = { 73, "eml_min_or_max",
  "C:/Program Files/MATLAB/R2013a/toolbox/eml/lib/matlab/eml/eml_min_or_max.m" };

static emlrtRSInfo i_emlrtRSI = { 88, "eml_min_or_max",
  "C:/Program Files/MATLAB/R2013a/toolbox/eml/lib/matlab/eml/eml_min_or_max.m" };

static emlrtRSInfo j_emlrtRSI = { 219, "eml_min_or_max",
  "C:/Program Files/MATLAB/R2013a/toolbox/eml/lib/matlab/eml/eml_min_or_max.m" };

static emlrtRSInfo k_emlrtRSI = { 192, "eml_min_or_max",
  "C:/Program Files/MATLAB/R2013a/toolbox/eml/lib/matlab/eml/eml_min_or_max.m" };

static emlrtRSInfo l_emlrtRSI = { 12, "eml_int_forloop_overflow_check",
  "C:/Program Files/MATLAB/R2013a/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"
};

static emlrtRSInfo m_emlrtRSI = { 51, "eml_int_forloop_overflow_check",
  "C:/Program Files/MATLAB/R2013a/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"
};

static emlrtRSInfo n_emlrtRSI = { 11, "eml_li_find",
  "C:/Program Files/MATLAB/R2013a/toolbox/eml/lib/matlab/eml/eml_li_find.m" };

static emlrtRSInfo o_emlrtRSI = { 14, "eml_li_find",
  "C:/Program Files/MATLAB/R2013a/toolbox/eml/lib/matlab/eml/eml_li_find.m" };

static emlrtRSInfo p_emlrtRSI = { 26, "eml_li_find",
  "C:/Program Files/MATLAB/R2013a/toolbox/eml/lib/matlab/eml/eml_li_find.m" };

static emlrtRSInfo q_emlrtRSI = { 39, "eml_li_find",
  "C:/Program Files/MATLAB/R2013a/toolbox/eml/lib/matlab/eml/eml_li_find.m" };

static emlrtRSInfo r_emlrtRSI = { 16, "max",
  "C:/Program Files/MATLAB/R2013a/toolbox/eml/lib/matlab/datafun/max.m" };

static emlrtMCInfo emlrtMCI = { 41, 9, "eml_min_or_max",
  "C:/Program Files/MATLAB/R2013a/toolbox/eml/lib/matlab/eml/eml_min_or_max.m" };

static emlrtMCInfo b_emlrtMCI = { 38, 19, "eml_min_or_max",
  "C:/Program Files/MATLAB/R2013a/toolbox/eml/lib/matlab/eml/eml_min_or_max.m" };

static emlrtMCInfo c_emlrtMCI = { 74, 9, "eml_min_or_max",
  "C:/Program Files/MATLAB/R2013a/toolbox/eml/lib/matlab/eml/eml_min_or_max.m" };

static emlrtMCInfo d_emlrtMCI = { 73, 19, "eml_min_or_max",
  "C:/Program Files/MATLAB/R2013a/toolbox/eml/lib/matlab/eml/eml_min_or_max.m" };

static emlrtMCInfo e_emlrtMCI = { 52, 9, "eml_int_forloop_overflow_check",
  "C:/Program Files/MATLAB/R2013a/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"
};

static emlrtMCInfo f_emlrtMCI = { 51, 15, "eml_int_forloop_overflow_check",
  "C:/Program Files/MATLAB/R2013a/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"
};

static emlrtMCInfo g_emlrtMCI = { 14, 5, "eml_li_find",
  "C:/Program Files/MATLAB/R2013a/toolbox/eml/lib/matlab/eml/eml_li_find.m" };

static emlrtRTEInfo b_emlrtRTEI = { 20, 9, "eml_li_find",
  "C:/Program Files/MATLAB/R2013a/toolbox/eml/lib/matlab/eml/eml_li_find.m" };

static emlrtRTEInfo e_emlrtRTEI = { 1, 26, "vibor_t",
  "C:/Users/LIK/Documents/MATLAB/Tube_dev/2be_0.660/Sim_Lib/vibor_t.m" };

static emlrtRTEInfo f_emlrtRTEI = { 6, 1, "vibor_t",
  "C:/Users/LIK/Documents/MATLAB/Tube_dev/2be_0.660/Sim_Lib/vibor_t.m" };

static emlrtRTEInfo g_emlrtRTEI = { 7, 1, "vibor_t",
  "C:/Users/LIK/Documents/MATLAB/Tube_dev/2be_0.660/Sim_Lib/vibor_t.m" };

static emlrtRTEInfo h_emlrtRTEI = { 8, 1, "vibor_t",
  "C:/Users/LIK/Documents/MATLAB/Tube_dev/2be_0.660/Sim_Lib/vibor_t.m" };

static emlrtRTEInfo i_emlrtRTEI = { 9, 1, "vibor_t",
  "C:/Users/LIK/Documents/MATLAB/Tube_dev/2be_0.660/Sim_Lib/vibor_t.m" };

static emlrtRTEInfo j_emlrtRTEI = { 11, 1, "vibor_t",
  "C:/Users/LIK/Documents/MATLAB/Tube_dev/2be_0.660/Sim_Lib/vibor_t.m" };

static emlrtRTEInfo k_emlrtRTEI = { 12, 1, "vibor_t",
  "C:/Users/LIK/Documents/MATLAB/Tube_dev/2be_0.660/Sim_Lib/vibor_t.m" };

static emlrtECInfo emlrtECI = { -1, 26, 46, "vibor_t",
  "C:/Users/LIK/Documents/MATLAB/Tube_dev/2be_0.660/Sim_Lib/vibor_t.m" };

static emlrtECInfo b_emlrtECI = { -1, 26, 21, "vibor_t",
  "C:/Users/LIK/Documents/MATLAB/Tube_dev/2be_0.660/Sim_Lib/vibor_t.m" };

static emlrtECInfo c_emlrtECI = { -1, 25, 34, "vibor_t",
  "C:/Users/LIK/Documents/MATLAB/Tube_dev/2be_0.660/Sim_Lib/vibor_t.m" };

static emlrtECInfo d_emlrtECI = { -1, 25, 21, "vibor_t",
  "C:/Users/LIK/Documents/MATLAB/Tube_dev/2be_0.660/Sim_Lib/vibor_t.m" };

static emlrtECInfo e_emlrtECI = { -1, 23, 9, "vibor_t",
  "C:/Users/LIK/Documents/MATLAB/Tube_dev/2be_0.660/Sim_Lib/vibor_t.m" };

static emlrtECInfo f_emlrtECI = { -1, 23, 29, "vibor_t",
  "C:/Users/LIK/Documents/MATLAB/Tube_dev/2be_0.660/Sim_Lib/vibor_t.m" };

static emlrtECInfo g_emlrtECI = { -1, 22, 9, "vibor_t",
  "C:/Users/LIK/Documents/MATLAB/Tube_dev/2be_0.660/Sim_Lib/vibor_t.m" };

static emlrtECInfo h_emlrtECI = { -1, 22, 29, "vibor_t",
  "C:/Users/LIK/Documents/MATLAB/Tube_dev/2be_0.660/Sim_Lib/vibor_t.m" };

static emlrtECInfo i_emlrtECI = { -1, 16, 14, "vibor_t",
  "C:/Users/LIK/Documents/MATLAB/Tube_dev/2be_0.660/Sim_Lib/vibor_t.m" };

static emlrtECInfo j_emlrtECI = { -1, 9, 6, "vibor_t",
  "C:/Users/LIK/Documents/MATLAB/Tube_dev/2be_0.660/Sim_Lib/vibor_t.m" };

static emlrtECInfo k_emlrtECI = { -1, 8, 6, "vibor_t",
  "C:/Users/LIK/Documents/MATLAB/Tube_dev/2be_0.660/Sim_Lib/vibor_t.m" };

static emlrtECInfo l_emlrtECI = { -1, 7, 5, "vibor_t",
  "C:/Users/LIK/Documents/MATLAB/Tube_dev/2be_0.660/Sim_Lib/vibor_t.m" };

static emlrtECInfo m_emlrtECI = { -1, 6, 5, "vibor_t",
  "C:/Users/LIK/Documents/MATLAB/Tube_dev/2be_0.660/Sim_Lib/vibor_t.m" };

static emlrtBCInfo emlrtBCI = { -1, -1, 6, 5, "Pc", "vibor_t",
  "C:/Users/LIK/Documents/MATLAB/Tube_dev/2be_0.660/Sim_Lib/vibor_t.m", 0 };

static emlrtDCInfo emlrtDCI = { 6, 5, "vibor_t",
  "C:/Users/LIK/Documents/MATLAB/Tube_dev/2be_0.660/Sim_Lib/vibor_t.m", 1 };

static emlrtBCInfo b_emlrtBCI = { -1, -1, 6, 15, "Pc", "vibor_t",
  "C:/Users/LIK/Documents/MATLAB/Tube_dev/2be_0.660/Sim_Lib/vibor_t.m", 0 };

static emlrtDCInfo b_emlrtDCI = { 6, 15, "vibor_t",
  "C:/Users/LIK/Documents/MATLAB/Tube_dev/2be_0.660/Sim_Lib/vibor_t.m", 1 };

static emlrtBCInfo c_emlrtBCI = { -1, -1, 7, 5, "Pg", "vibor_t",
  "C:/Users/LIK/Documents/MATLAB/Tube_dev/2be_0.660/Sim_Lib/vibor_t.m", 0 };

static emlrtDCInfo c_emlrtDCI = { 7, 5, "vibor_t",
  "C:/Users/LIK/Documents/MATLAB/Tube_dev/2be_0.660/Sim_Lib/vibor_t.m", 1 };

static emlrtBCInfo d_emlrtBCI = { -1, -1, 7, 15, "Pg", "vibor_t",
  "C:/Users/LIK/Documents/MATLAB/Tube_dev/2be_0.660/Sim_Lib/vibor_t.m", 0 };

static emlrtDCInfo d_emlrtDCI = { 7, 15, "vibor_t",
  "C:/Users/LIK/Documents/MATLAB/Tube_dev/2be_0.660/Sim_Lib/vibor_t.m", 1 };

static emlrtBCInfo e_emlrtBCI = { -1, -1, 8, 6, "Pc", "vibor_t",
  "C:/Users/LIK/Documents/MATLAB/Tube_dev/2be_0.660/Sim_Lib/vibor_t.m", 0 };

static emlrtDCInfo e_emlrtDCI = { 8, 6, "vibor_t",
  "C:/Users/LIK/Documents/MATLAB/Tube_dev/2be_0.660/Sim_Lib/vibor_t.m", 1 };

static emlrtBCInfo f_emlrtBCI = { -1, -1, 8, 19, "Pw", "vibor_t",
  "C:/Users/LIK/Documents/MATLAB/Tube_dev/2be_0.660/Sim_Lib/vibor_t.m", 0 };

static emlrtDCInfo f_emlrtDCI = { 8, 19, "vibor_t",
  "C:/Users/LIK/Documents/MATLAB/Tube_dev/2be_0.660/Sim_Lib/vibor_t.m", 1 };

static emlrtBCInfo g_emlrtBCI = { -1, -1, 9, 6, "Pg", "vibor_t",
  "C:/Users/LIK/Documents/MATLAB/Tube_dev/2be_0.660/Sim_Lib/vibor_t.m", 0 };

static emlrtDCInfo g_emlrtDCI = { 9, 6, "vibor_t",
  "C:/Users/LIK/Documents/MATLAB/Tube_dev/2be_0.660/Sim_Lib/vibor_t.m", 1 };

static emlrtBCInfo h_emlrtBCI = { -1, -1, 9, 19, "Pw", "vibor_t",
  "C:/Users/LIK/Documents/MATLAB/Tube_dev/2be_0.660/Sim_Lib/vibor_t.m", 0 };

static emlrtDCInfo h_emlrtDCI = { 9, 19, "vibor_t",
  "C:/Users/LIK/Documents/MATLAB/Tube_dev/2be_0.660/Sim_Lib/vibor_t.m", 1 };

static emlrtBCInfo i_emlrtBCI = { -1, -1, 22, 9, "dVc", "vibor_t",
  "C:/Users/LIK/Documents/MATLAB/Tube_dev/2be_0.660/Sim_Lib/vibor_t.m", 0 };

static emlrtDCInfo i_emlrtDCI = { 22, 9, "vibor_t",
  "C:/Users/LIK/Documents/MATLAB/Tube_dev/2be_0.660/Sim_Lib/vibor_t.m", 1 };

static emlrtBCInfo j_emlrtBCI = { -1, -1, 22, 29, "dVc", "vibor_t",
  "C:/Users/LIK/Documents/MATLAB/Tube_dev/2be_0.660/Sim_Lib/vibor_t.m", 0 };

static emlrtDCInfo j_emlrtDCI = { 22, 29, "vibor_t",
  "C:/Users/LIK/Documents/MATLAB/Tube_dev/2be_0.660/Sim_Lib/vibor_t.m", 1 };

static emlrtBCInfo k_emlrtBCI = { -1, -1, 23, 9, "dVg", "vibor_t",
  "C:/Users/LIK/Documents/MATLAB/Tube_dev/2be_0.660/Sim_Lib/vibor_t.m", 0 };

static emlrtDCInfo k_emlrtDCI = { 23, 9, "vibor_t",
  "C:/Users/LIK/Documents/MATLAB/Tube_dev/2be_0.660/Sim_Lib/vibor_t.m", 1 };

static emlrtBCInfo l_emlrtBCI = { -1, -1, 23, 29, "dVg", "vibor_t",
  "C:/Users/LIK/Documents/MATLAB/Tube_dev/2be_0.660/Sim_Lib/vibor_t.m", 0 };

static emlrtDCInfo l_emlrtDCI = { 23, 29, "vibor_t",
  "C:/Users/LIK/Documents/MATLAB/Tube_dev/2be_0.660/Sim_Lib/vibor_t.m", 1 };

static emlrtBCInfo m_emlrtBCI = { -1, -1, 26, 57, "dVg", "vibor_t",
  "C:/Users/LIK/Documents/MATLAB/Tube_dev/2be_0.660/Sim_Lib/vibor_t.m", 0 };

static emlrtDCInfo m_emlrtDCI = { 26, 57, "vibor_t",
  "C:/Users/LIK/Documents/MATLAB/Tube_dev/2be_0.660/Sim_Lib/vibor_t.m", 1 };

static emlrtBCInfo n_emlrtBCI = { -1, -1, 26, 32, "dVc", "vibor_t",
  "C:/Users/LIK/Documents/MATLAB/Tube_dev/2be_0.660/Sim_Lib/vibor_t.m", 0 };

static emlrtDCInfo n_emlrtDCI = { 26, 32, "vibor_t",
  "C:/Users/LIK/Documents/MATLAB/Tube_dev/2be_0.660/Sim_Lib/vibor_t.m", 1 };

static emlrtBCInfo o_emlrtBCI = { -1, -1, 18, 12, "MdS", "vibor_t",
  "C:/Users/LIK/Documents/MATLAB/Tube_dev/2be_0.660/Sim_Lib/vibor_t.m", 0 };

/* Function Declarations */
static const mxArray *b_message(const mxArray *b, const mxArray *c, emlrtMCInfo *
  location);
static void check_forloop_overflow_error(void);
static void eml_li_find(const emxArray_boolean_T *x, emxArray_int32_T *y);

/* Function Definitions */
static const mxArray *b_message(const mxArray *b, const mxArray *c, emlrtMCInfo *
  location)
{
  const mxArray *pArrays[2];
  const mxArray *m7;
  pArrays[0] = b;
  pArrays[1] = c;
  return emlrtCallMATLABR2012b(emlrtRootTLSGlobal, 1, &m7, 2, pArrays, "message",
    TRUE, location);
}

static void check_forloop_overflow_error(void)
{
  const mxArray *y;
  static const int32_T iv0[2] = { 1, 34 };

  const mxArray *m0;
  char_T cv0[34];
  int32_T i;
  static const char_T cv1[34] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o', 'o',
    'l', 'b', 'o', 'x', ':', 'i', 'n', 't', '_', 'f', 'o', 'r', 'l', 'o', 'o',
    'p', '_', 'o', 'v', 'e', 'r', 'f', 'l', 'o', 'w' };

  const mxArray *b_y;
  static const int32_T iv1[2] = { 1, 23 };

  char_T cv2[23];
  static const char_T cv3[23] = { 'c', 'o', 'd', 'e', 'r', '.', 'i', 'n', 't',
    'e', 'r', 'n', 'a', 'l', '.', 'i', 'n', 'd', 'e', 'x', 'I', 'n', 't' };

  emlrtPushRtStackR2012b(&m_emlrtRSI, emlrtRootTLSGlobal);
  y = NULL;
  m0 = mxCreateCharArray(2, iv0);
  for (i = 0; i < 34; i++) {
    cv0[i] = cv1[i];
  }

  emlrtInitCharArrayR2013a(emlrtRootTLSGlobal, 34, m0, cv0);
  emlrtAssign(&y, m0);
  b_y = NULL;
  m0 = mxCreateCharArray(2, iv1);
  for (i = 0; i < 23; i++) {
    cv2[i] = cv3[i];
  }

  emlrtInitCharArrayR2013a(emlrtRootTLSGlobal, 23, m0, cv2);
  emlrtAssign(&b_y, m0);
  error(b_message(y, b_y, &e_emlrtMCI), &f_emlrtMCI);
  emlrtPopRtStackR2012b(&m_emlrtRSI, emlrtRootTLSGlobal);
}

static void eml_li_find(const emxArray_boolean_T *x, emxArray_int32_T *y)
{
  int32_T k;
  boolean_T overflow;
  int32_T i;
  const mxArray *b_y;
  const mxArray *m1;
  int32_T j;
  emlrtPushRtStackR2012b(&n_emlrtRSI, emlrtRootTLSGlobal);
  k = 0;
  emlrtPushRtStackR2012b(&q_emlrtRSI, emlrtRootTLSGlobal);
  if (1 > x->size[0]) {
    overflow = FALSE;
  } else {
    overflow = (x->size[0] > 2147483646);
  }

  if (overflow) {
    emlrtPushRtStackR2012b(&l_emlrtRSI, emlrtRootTLSGlobal);
    check_forloop_overflow_error();
    emlrtPopRtStackR2012b(&l_emlrtRSI, emlrtRootTLSGlobal);
  }

  emlrtPopRtStackR2012b(&q_emlrtRSI, emlrtRootTLSGlobal);
  for (i = 1; i <= x->size[0]; i++) {
    if (x->data[i - 1]) {
      k++;
    }
  }

  emlrtPopRtStackR2012b(&n_emlrtRSI, emlrtRootTLSGlobal);
  if (k <= x->size[0]) {
  } else {
    emlrtPushRtStackR2012b(&o_emlrtRSI, emlrtRootTLSGlobal);
    b_y = NULL;
    m1 = mxCreateString("Assertion failed.");
    emlrtAssign(&b_y, m1);
    error(b_y, &g_emlrtMCI);
    emlrtPopRtStackR2012b(&o_emlrtRSI, emlrtRootTLSGlobal);
  }

  j = y->size[0];
  y->size[0] = k;
  emxEnsureCapacity((emxArray__common *)y, j, (int32_T)sizeof(int32_T),
                    &b_emlrtRTEI);
  j = 0;
  emlrtPushRtStackR2012b(&p_emlrtRSI, emlrtRootTLSGlobal);
  if (1 > x->size[0]) {
    overflow = FALSE;
  } else {
    overflow = (x->size[0] > 2147483646);
  }

  if (overflow) {
    emlrtPushRtStackR2012b(&l_emlrtRSI, emlrtRootTLSGlobal);
    check_forloop_overflow_error();
    emlrtPopRtStackR2012b(&l_emlrtRSI, emlrtRootTLSGlobal);
  }

  emlrtPopRtStackR2012b(&p_emlrtRSI, emlrtRootTLSGlobal);
  for (i = 1; i <= x->size[0]; i++) {
    if (x->data[i - 1]) {
      y->data[j] = i;
      j++;
    }
  }
}

boolean_T vibor_t(real_T *ndt, real_T fl, const emxArray_real_T *Pc, const
                  emxArray_real_T *Pg, const emxArray_real_T *Pw, const
                  b_struct_T *PR, const struct_T *RC, real_T dt, const
                  emxArray_real_T *Sw1, const emxArray_real_T *Sw2, const
                  emxArray_real_T *CL, const emxArray_real_T *GL, const
                  emxArray_real_T *dVc, const emxArray_real_T *dVg, const
                  emxArray_real_T *W1C, const emxArray_real_T *W1G, const
                  emxArray_real_T *WoC, const emxArray_real_T *WoG, real_T i,
                  real_T *j_ndt)
{
  boolean_T fl1;
  int32_T b_i;
  int32_T ix;
  int32_T i6;
  real_T j_old;
  int32_T i7;
  emxArray_real_T *dPc;
  emxArray_real_T *dPg;
  emxArray_real_T *dPwc;
  emxArray_real_T *dPwg;
  emxArray_boolean_T *v1;
  emxArray_boolean_T *v2;
  boolean_T overflow;
  const mxArray *y;
  static const int32_T iv15[2] = { 1, 36 };

  const mxArray *m8;
  char_T cv6[36];
  static const char_T cv7[36] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o', 'o',
    'l', 'b', 'o', 'x', ':', 'a', 'u', 't', 'o', 'D', 'i', 'm', 'I', 'n', 'c',
    'o', 'm', 'p', 'a', 't', 'i', 'b', 'i', 'l', 'i', 't', 'y' };

  const mxArray *b_y;
  static const int32_T iv16[2] = { 1, 39 };

  char_T cv8[39];
  static const char_T cv9[39] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o', 'o',
    'l', 'b', 'o', 'x', ':', 'e', 'm', 'l', '_', 'm', 'i', 'n', '_', 'o', 'r',
    '_', 'm', 'a', 'x', '_', 'v', 'a', 'r', 'D', 'i', 'm', 'Z', 'e', 'r', 'o' };

  real_T mtmp;
  boolean_T exitg5;
  emxArray_boolean_T *r1;
  emxArray_real_T *r2;
  emxArray_real_T *r3;
  emxArray_real_T *r4;
  emxArray_real_T *r5;
  emxArray_int32_T *r6;
  emxArray_real_T *r7;
  emxArray_real_T *b_W1C;
  emxArray_real_T *b_dVc;
  emxArray_real_T *b_W1G;
  emxArray_real_T *b_dVg;
  emxArray_real_T *r8;
  emxArray_real_T *b_CL;
  emxArray_real_T *r9;
  emxArray_real_T *b_GL;
  emxArray_real_T *r10;
  emxArray_boolean_T *r11;
  emxArray_real_T *b_Sw2;
  const mxArray *c_y;
  static const int32_T iv17[2] = { 1, 36 };

  const mxArray *d_y;
  static const int32_T iv18[2] = { 1, 39 };

  boolean_T exitg4;
  const mxArray *e_y;
  static const int32_T iv19[2] = { 1, 36 };

  const mxArray *f_y;
  static const int32_T iv20[2] = { 1, 39 };

  boolean_T exitg3;
  const mxArray *g_y;
  static const int32_T iv21[2] = { 1, 36 };

  const mxArray *h_y;
  static const int32_T iv22[2] = { 1, 39 };

  boolean_T exitg2;
  boolean_T exitg1;
  emlrtHeapReferenceStackEnterFcnR2012b(emlrtRootTLSGlobal);
  b_i = RC->Cc->size[0];
  for (ix = 0; ix < b_i; ix++) {
    i6 = Pc->size[0];
    j_old = RC->Cc->data[ix];
    i7 = (int32_T)emlrtIntegerCheckFastR2012b(j_old, &emlrtDCI,
      emlrtRootTLSGlobal);
    emlrtDynamicBoundsCheckFastR2012b(i7, 1, i6, &emlrtBCI, emlrtRootTLSGlobal);
  }

  b_i = RC->Cr->size[0];
  for (ix = 0; ix < b_i; ix++) {
    i6 = Pc->size[0];
    j_old = RC->Cr->data[ix];
    i7 = (int32_T)emlrtIntegerCheckFastR2012b(j_old, &b_emlrtDCI,
      emlrtRootTLSGlobal);
    emlrtDynamicBoundsCheckFastR2012b(i7, 1, i6, &b_emlrtBCI, emlrtRootTLSGlobal);
  }

  emxInit_real_T(&dPc, 1, &f_emlrtRTEI, TRUE);
  ix = RC->Cc->size[0];
  i6 = RC->Cr->size[0];
  emlrtSizeEqCheck1DFastR2012b(ix, i6, &m_emlrtECI, emlrtRootTLSGlobal);
  ix = dPc->size[0];
  dPc->size[0] = RC->Cc->size[0];
  emxEnsureCapacity((emxArray__common *)dPc, ix, (int32_T)sizeof(real_T),
                    &e_emlrtRTEI);
  b_i = RC->Cc->size[0];
  for (ix = 0; ix < b_i; ix++) {
    dPc->data[ix] = Pc->data[(int32_T)RC->Cc->data[ix] - 1] - Pc->data[(int32_T)
      RC->Cr->data[ix] - 1];
  }

  b_i = RC->Gc->size[0];
  for (ix = 0; ix < b_i; ix++) {
    i6 = Pg->size[0];
    j_old = RC->Gc->data[ix];
    i7 = (int32_T)emlrtIntegerCheckFastR2012b(j_old, &c_emlrtDCI,
      emlrtRootTLSGlobal);
    emlrtDynamicBoundsCheckFastR2012b(i7, 1, i6, &c_emlrtBCI, emlrtRootTLSGlobal);
  }

  b_i = RC->Gr->size[0];
  for (ix = 0; ix < b_i; ix++) {
    i6 = Pg->size[0];
    j_old = RC->Gr->data[ix];
    i7 = (int32_T)emlrtIntegerCheckFastR2012b(j_old, &d_emlrtDCI,
      emlrtRootTLSGlobal);
    emlrtDynamicBoundsCheckFastR2012b(i7, 1, i6, &d_emlrtBCI, emlrtRootTLSGlobal);
  }

  emxInit_real_T(&dPg, 1, &g_emlrtRTEI, TRUE);
  ix = RC->Gc->size[0];
  i6 = RC->Gr->size[0];
  emlrtSizeEqCheck1DFastR2012b(ix, i6, &l_emlrtECI, emlrtRootTLSGlobal);
  ix = dPg->size[0];
  dPg->size[0] = RC->Gc->size[0];
  emxEnsureCapacity((emxArray__common *)dPg, ix, (int32_T)sizeof(real_T),
                    &e_emlrtRTEI);
  b_i = RC->Gc->size[0];
  for (ix = 0; ix < b_i; ix++) {
    dPg->data[ix] = Pg->data[(int32_T)RC->Gc->data[ix] - 1] - Pg->data[(int32_T)
      RC->Gr->data[ix] - 1];
  }

  b_i = WoC->size[0] - 1;
  for (ix = 0; ix <= b_i; ix++) {
    i6 = Pc->size[0];
    j_old = WoC->data[ix];
    i7 = (int32_T)emlrtIntegerCheckFastR2012b(j_old, &e_emlrtDCI,
      emlrtRootTLSGlobal);
    emlrtDynamicBoundsCheckFastR2012b(i7, 1, i6, &e_emlrtBCI, emlrtRootTLSGlobal);
  }

  b_i = WoC->size[0] - 1;
  for (ix = 0; ix <= b_i; ix++) {
    i6 = Pw->size[0];
    j_old = WoC->data[ix + (WoC->size[0] << 1)];
    i7 = (int32_T)emlrtIntegerCheckFastR2012b(j_old, &f_emlrtDCI,
      emlrtRootTLSGlobal);
    emlrtDynamicBoundsCheckFastR2012b(i7, 1, i6, &f_emlrtBCI, emlrtRootTLSGlobal);
  }

  emxInit_real_T(&dPwc, 1, &h_emlrtRTEI, TRUE);
  ix = WoC->size[0];
  i6 = WoC->size[0];
  emlrtSizeEqCheck1DFastR2012b(ix, i6, &k_emlrtECI, emlrtRootTLSGlobal);
  b_i = WoC->size[0];
  ix = dPwc->size[0];
  dPwc->size[0] = b_i;
  emxEnsureCapacity((emxArray__common *)dPwc, ix, (int32_T)sizeof(real_T),
                    &e_emlrtRTEI);
  for (ix = 0; ix < b_i; ix++) {
    dPwc->data[ix] = Pc->data[(int32_T)WoC->data[ix] - 1] - Pw->data[(int32_T)
      WoC->data[ix + (WoC->size[0] << 1)] - 1];
  }

  b_i = WoG->size[0] - 1;
  for (ix = 0; ix <= b_i; ix++) {
    i6 = Pg->size[0];
    j_old = WoG->data[ix];
    i7 = (int32_T)emlrtIntegerCheckFastR2012b(j_old, &g_emlrtDCI,
      emlrtRootTLSGlobal);
    emlrtDynamicBoundsCheckFastR2012b(i7, 1, i6, &g_emlrtBCI, emlrtRootTLSGlobal);
  }

  b_i = WoG->size[0] - 1;
  for (ix = 0; ix <= b_i; ix++) {
    i6 = Pw->size[0];
    j_old = WoG->data[ix + (WoG->size[0] << 1)];
    i7 = (int32_T)emlrtIntegerCheckFastR2012b(j_old, &h_emlrtDCI,
      emlrtRootTLSGlobal);
    emlrtDynamicBoundsCheckFastR2012b(i7, 1, i6, &h_emlrtBCI, emlrtRootTLSGlobal);
  }

  emxInit_real_T(&dPwg, 1, &i_emlrtRTEI, TRUE);
  ix = WoG->size[0];
  i6 = WoG->size[0];
  emlrtSizeEqCheck1DFastR2012b(ix, i6, &j_emlrtECI, emlrtRootTLSGlobal);
  b_i = WoG->size[0];
  ix = dPwg->size[0];
  dPwg->size[0] = b_i;
  emxEnsureCapacity((emxArray__common *)dPwg, ix, (int32_T)sizeof(real_T),
                    &e_emlrtRTEI);
  for (ix = 0; ix < b_i; ix++) {
    dPwg->data[ix] = Pg->data[(int32_T)WoG->data[ix] - 1] - Pw->data[(int32_T)
      WoG->data[ix + (WoG->size[0] << 1)] - 1];
  }

  emxInit_boolean_T(&v1, 1, &j_emlrtRTEI, TRUE);
  ix = v1->size[0];
  v1->size[0] = dPc->size[0];
  emxEnsureCapacity((emxArray__common *)v1, ix, (int32_T)sizeof(boolean_T),
                    &e_emlrtRTEI);
  b_i = dPc->size[0];
  for (ix = 0; ix < b_i; ix++) {
    v1->data[ix] = (dPc->data[ix] > 0.0);
  }

  emxInit_boolean_T(&v2, 1, &k_emlrtRTEI, TRUE);
  ix = v2->size[0];
  v2->size[0] = dPg->size[0];
  emxEnsureCapacity((emxArray__common *)v2, ix, (int32_T)sizeof(boolean_T),
                    &e_emlrtRTEI);
  b_i = dPg->size[0];
  for (ix = 0; ix < b_i; ix++) {
    v2->data[ix] = (dPg->data[ix] > 0.0);
  }

  if (fl == 0.0) {
    emlrtPushRtStackR2012b(&emlrtRSI, emlrtRootTLSGlobal);
    emlrtPushRtStackR2012b(&e_emlrtRSI, emlrtRootTLSGlobal);
    emlrtPushRtStackR2012b(&f_emlrtRSI, emlrtRootTLSGlobal);
    if ((Sw1->size[0] == 1) || (Sw1->size[0] != 1)) {
      overflow = TRUE;
    } else {
      overflow = FALSE;
    }

    if (overflow) {
    } else {
      emlrtPushRtStackR2012b(&g_emlrtRSI, emlrtRootTLSGlobal);
      y = NULL;
      m8 = mxCreateCharArray(2, iv15);
      for (b_i = 0; b_i < 36; b_i++) {
        cv6[b_i] = cv7[b_i];
      }

      emlrtInitCharArrayR2013a(emlrtRootTLSGlobal, 36, m8, cv6);
      emlrtAssign(&y, m8);
      error(message(y, &emlrtMCI), &b_emlrtMCI);
      emlrtPopRtStackR2012b(&g_emlrtRSI, emlrtRootTLSGlobal);
    }

    if (Sw1->size[0] > 0) {
    } else {
      emlrtPushRtStackR2012b(&h_emlrtRSI, emlrtRootTLSGlobal);
      b_y = NULL;
      m8 = mxCreateCharArray(2, iv16);
      for (b_i = 0; b_i < 39; b_i++) {
        cv8[b_i] = cv9[b_i];
      }

      emlrtInitCharArrayR2013a(emlrtRootTLSGlobal, 39, m8, cv8);
      emlrtAssign(&b_y, m8);
      error(message(b_y, &c_emlrtMCI), &d_emlrtMCI);
      emlrtPopRtStackR2012b(&h_emlrtRSI, emlrtRootTLSGlobal);
    }

    emlrtPushRtStackR2012b(&i_emlrtRSI, emlrtRootTLSGlobal);
    b_i = 1;
    mtmp = Sw1->data[0];
    if (Sw1->size[0] > 1) {
      if (muDoubleScalarIsNaN(Sw1->data[0])) {
        emlrtPushRtStackR2012b(&k_emlrtRSI, emlrtRootTLSGlobal);
        if (2 > Sw1->size[0]) {
          overflow = FALSE;
        } else {
          overflow = (Sw1->size[0] > 2147483646);
        }

        if (overflow) {
          emlrtPushRtStackR2012b(&l_emlrtRSI, emlrtRootTLSGlobal);
          check_forloop_overflow_error();
          emlrtPopRtStackR2012b(&l_emlrtRSI, emlrtRootTLSGlobal);
        }

        emlrtPopRtStackR2012b(&k_emlrtRSI, emlrtRootTLSGlobal);
        ix = 2;
        exitg5 = FALSE;
        while ((exitg5 == FALSE) && (ix <= Sw1->size[0])) {
          b_i = ix;
          if (!muDoubleScalarIsNaN(Sw1->data[ix - 1])) {
            mtmp = Sw1->data[ix - 1];
            exitg5 = TRUE;
          } else {
            ix++;
          }
        }
      }

      if (b_i < Sw1->size[0]) {
        emlrtPushRtStackR2012b(&j_emlrtRSI, emlrtRootTLSGlobal);
        if (b_i + 1 > Sw1->size[0]) {
          overflow = FALSE;
        } else {
          overflow = (Sw1->size[0] > 2147483646);
        }

        if (overflow) {
          emlrtPushRtStackR2012b(&l_emlrtRSI, emlrtRootTLSGlobal);
          check_forloop_overflow_error();
          emlrtPopRtStackR2012b(&l_emlrtRSI, emlrtRootTLSGlobal);
        }

        emlrtPopRtStackR2012b(&j_emlrtRSI, emlrtRootTLSGlobal);
        while (b_i + 1 <= Sw1->size[0]) {
          if (Sw1->data[b_i] < mtmp) {
            mtmp = Sw1->data[b_i];
          }

          b_i++;
        }
      }
    }

    emlrtPopRtStackR2012b(&i_emlrtRSI, emlrtRootTLSGlobal);
    emlrtPopRtStackR2012b(&f_emlrtRSI, emlrtRootTLSGlobal);
    emlrtPopRtStackR2012b(&e_emlrtRSI, emlrtRootTLSGlobal);
    emlrtPopRtStackR2012b(&emlrtRSI, emlrtRootTLSGlobal);
    emxInit_boolean_T(&r1, 1, &e_emlrtRTEI, TRUE);
    emxInit_real_T(&r2, 1, &e_emlrtRTEI, TRUE);
    emxInit_real_T(&r3, 1, &e_emlrtRTEI, TRUE);
    emxInit_real_T(&r4, 1, &e_emlrtRTEI, TRUE);
    emxInit_real_T(&r5, 1, &e_emlrtRTEI, TRUE);
    emxInit_int32_T(&r6, 1, &e_emlrtRTEI, TRUE);
    emxInit_real_T(&r7, 1, &e_emlrtRTEI, TRUE);
    emxInit_real_T(&b_W1C, 1, &e_emlrtRTEI, TRUE);
    emxInit_real_T(&b_dVc, 1, &e_emlrtRTEI, TRUE);
    emxInit_real_T(&b_W1G, 1, &e_emlrtRTEI, TRUE);
    emxInit_real_T(&b_dVg, 1, &e_emlrtRTEI, TRUE);
    emxInit_real_T(&r8, 1, &e_emlrtRTEI, TRUE);
    emxInit_real_T(&b_CL, 1, &e_emlrtRTEI, TRUE);
    emxInit_real_T(&r9, 1, &e_emlrtRTEI, TRUE);
    emxInit_real_T(&b_GL, 1, &e_emlrtRTEI, TRUE);
    emxInit_real_T(&r10, 1, &e_emlrtRTEI, TRUE);
    emxInit_boolean_T(&r11, 1, &e_emlrtRTEI, TRUE);
    emxInit_real_T(&b_Sw2, 1, &e_emlrtRTEI, TRUE);
    if ((mtmp > PR->Sc2) && (i != 1.0)) {
      ix = Sw2->size[0];
      i6 = Sw1->size[0];
      emlrtSizeEqCheck1DFastR2012b(ix, i6, &i_emlrtECI, emlrtRootTLSGlobal);
      ix = b_Sw2->size[0];
      b_Sw2->size[0] = Sw2->size[0];
      emxEnsureCapacity((emxArray__common *)b_Sw2, ix, (int32_T)sizeof(real_T),
                        &e_emlrtRTEI);
      b_i = Sw2->size[0];
      for (ix = 0; ix < b_i; ix++) {
        b_Sw2->data[ix] = Sw2->data[ix] - Sw1->data[ix];
      }

      b_abs(b_Sw2, dPc);
      emlrtPushRtStackR2012b(&b_emlrtRSI, emlrtRootTLSGlobal);
      ix = r1->size[0];
      r1->size[0] = dPc->size[0];
      emxEnsureCapacity((emxArray__common *)r1, ix, (int32_T)sizeof(boolean_T),
                        &e_emlrtRTEI);
      b_i = dPc->size[0];
      for (ix = 0; ix < b_i; ix++) {
        r1->data[ix] = muDoubleScalarIsNaN(dPc->data[ix]);
      }

      ix = r11->size[0];
      r11->size[0] = r1->size[0];
      emxEnsureCapacity((emxArray__common *)r11, ix, (int32_T)sizeof(boolean_T),
                        &e_emlrtRTEI);
      b_i = r1->size[0];
      for (ix = 0; ix < b_i; ix++) {
        r11->data[ix] = ((r1->data[ix] == 0) == 1);
      }

      eml_li_find(r11, r6);
      ix = dPg->size[0];
      dPg->size[0] = r6->size[0];
      emxEnsureCapacity((emxArray__common *)dPg, ix, (int32_T)sizeof(real_T),
                        &e_emlrtRTEI);
      b_i = r6->size[0];
      for (ix = 0; ix < b_i; ix++) {
        i6 = dPc->size[0];
        i7 = r6->data[ix];
        dPg->data[ix] = dPc->data[emlrtDynamicBoundsCheckFastR2012b(i7, 1, i6,
          &o_emlrtBCI, emlrtRootTLSGlobal) - 1];
      }

      emlrtPushRtStackR2012b(&r_emlrtRSI, emlrtRootTLSGlobal);
      emlrtPushRtStackR2012b(&f_emlrtRSI, emlrtRootTLSGlobal);
      if ((dPg->size[0] == 1) || (dPg->size[0] != 1)) {
        overflow = TRUE;
      } else {
        overflow = FALSE;
      }

      if (overflow) {
      } else {
        emlrtPushRtStackR2012b(&g_emlrtRSI, emlrtRootTLSGlobal);
        c_y = NULL;
        m8 = mxCreateCharArray(2, iv17);
        for (b_i = 0; b_i < 36; b_i++) {
          cv6[b_i] = cv7[b_i];
        }

        emlrtInitCharArrayR2013a(emlrtRootTLSGlobal, 36, m8, cv6);
        emlrtAssign(&c_y, m8);
        error(message(c_y, &emlrtMCI), &b_emlrtMCI);
        emlrtPopRtStackR2012b(&g_emlrtRSI, emlrtRootTLSGlobal);
      }

      if (dPg->size[0] > 0) {
      } else {
        emlrtPushRtStackR2012b(&h_emlrtRSI, emlrtRootTLSGlobal);
        d_y = NULL;
        m8 = mxCreateCharArray(2, iv18);
        for (b_i = 0; b_i < 39; b_i++) {
          cv8[b_i] = cv9[b_i];
        }

        emlrtInitCharArrayR2013a(emlrtRootTLSGlobal, 39, m8, cv8);
        emlrtAssign(&d_y, m8);
        error(message(d_y, &c_emlrtMCI), &d_emlrtMCI);
        emlrtPopRtStackR2012b(&h_emlrtRSI, emlrtRootTLSGlobal);
      }

      emlrtPushRtStackR2012b(&i_emlrtRSI, emlrtRootTLSGlobal);
      b_i = 1;
      mtmp = dPg->data[0];
      if (dPg->size[0] > 1) {
        if (muDoubleScalarIsNaN(dPg->data[0])) {
          emlrtPushRtStackR2012b(&k_emlrtRSI, emlrtRootTLSGlobal);
          if (2 > dPg->size[0]) {
            overflow = FALSE;
          } else {
            overflow = (dPg->size[0] > 2147483646);
          }

          if (overflow) {
            emlrtPushRtStackR2012b(&l_emlrtRSI, emlrtRootTLSGlobal);
            check_forloop_overflow_error();
            emlrtPopRtStackR2012b(&l_emlrtRSI, emlrtRootTLSGlobal);
          }

          emlrtPopRtStackR2012b(&k_emlrtRSI, emlrtRootTLSGlobal);
          ix = 2;
          exitg4 = FALSE;
          while ((exitg4 == FALSE) && (ix <= dPg->size[0])) {
            b_i = ix;
            if (!muDoubleScalarIsNaN(dPg->data[ix - 1])) {
              mtmp = dPg->data[ix - 1];
              exitg4 = TRUE;
            } else {
              ix++;
            }
          }
        }

        if (b_i < dPg->size[0]) {
          emlrtPushRtStackR2012b(&j_emlrtRSI, emlrtRootTLSGlobal);
          if (b_i + 1 > dPg->size[0]) {
            overflow = FALSE;
          } else {
            overflow = (dPg->size[0] > 2147483646);
          }

          if (overflow) {
            emlrtPushRtStackR2012b(&l_emlrtRSI, emlrtRootTLSGlobal);
            check_forloop_overflow_error();
            emlrtPopRtStackR2012b(&l_emlrtRSI, emlrtRootTLSGlobal);
          }

          emlrtPopRtStackR2012b(&j_emlrtRSI, emlrtRootTLSGlobal);
          while (b_i + 1 <= dPg->size[0]) {
            if (dPg->data[b_i] > mtmp) {
              mtmp = dPg->data[b_i];
            }

            b_i++;
          }
        }
      }

      emlrtPopRtStackR2012b(&i_emlrtRSI, emlrtRootTLSGlobal);
      emlrtPopRtStackR2012b(&f_emlrtRSI, emlrtRootTLSGlobal);
      emlrtPopRtStackR2012b(&r_emlrtRSI, emlrtRootTLSGlobal);
      emlrtPopRtStackR2012b(&b_emlrtRSI, emlrtRootTLSGlobal);
      *ndt = muDoubleScalarCeil(*ndt * (mtmp / 0.02));
    } else {
      b_i = RC->Cc->size[0];
      for (ix = 0; ix < b_i; ix++) {
        i6 = dVc->size[0];
        j_old = RC->Cc->data[ix];
        i7 = (int32_T)emlrtIntegerCheckFastR2012b(j_old, &i_emlrtDCI,
          emlrtRootTLSGlobal);
        emlrtDynamicBoundsCheckFastR2012b(i7, 1, i6, &i_emlrtBCI,
          emlrtRootTLSGlobal);
      }

      ix = r1->size[0];
      r1->size[0] = v1->size[0];
      emxEnsureCapacity((emxArray__common *)r1, ix, (int32_T)sizeof(boolean_T),
                        &e_emlrtRTEI);
      b_i = v1->size[0];
      for (ix = 0; ix < b_i; ix++) {
        r1->data[ix] = (v1->data[ix] == 0);
      }

      ix = RC->Cc->size[0];
      i6 = r1->size[0];
      emlrtSizeEqCheck1DFastR2012b(ix, i6, &g_emlrtECI, emlrtRootTLSGlobal);
      b_i = RC->Cr->size[0];
      for (ix = 0; ix < b_i; ix++) {
        i6 = dVc->size[0];
        j_old = RC->Cr->data[ix];
        i7 = (int32_T)emlrtIntegerCheckFastR2012b(j_old, &j_emlrtDCI,
          emlrtRootTLSGlobal);
        emlrtDynamicBoundsCheckFastR2012b(i7, 1, i6, &j_emlrtBCI,
          emlrtRootTLSGlobal);
      }

      ix = RC->Cr->size[0];
      i6 = v1->size[0];
      emlrtSizeEqCheck1DFastR2012b(ix, i6, &h_emlrtECI, emlrtRootTLSGlobal);
      ix = r2->size[0];
      r2->size[0] = RC->Cc->size[0];
      emxEnsureCapacity((emxArray__common *)r2, ix, (int32_T)sizeof(real_T),
                        &e_emlrtRTEI);
      b_i = RC->Cc->size[0];
      for (ix = 0; ix < b_i; ix++) {
        r2->data[ix] = dVc->data[(int32_T)RC->Cc->data[ix] - 1] * (real_T)
          r1->data[ix];
      }

      ix = r3->size[0];
      r3->size[0] = RC->Cr->size[0];
      emxEnsureCapacity((emxArray__common *)r3, ix, (int32_T)sizeof(real_T),
                        &e_emlrtRTEI);
      b_i = RC->Cr->size[0];
      for (ix = 0; ix < b_i; ix++) {
        r3->data[ix] = dVc->data[(int32_T)RC->Cr->data[ix] - 1] * (real_T)
          v1->data[ix];
      }

      ix = r2->size[0];
      i6 = r3->size[0];
      emlrtSizeEqCheck1DFastR2012b(ix, i6, &g_emlrtECI, emlrtRootTLSGlobal);
      b_i = RC->Gc->size[0];
      for (ix = 0; ix < b_i; ix++) {
        i6 = dVg->size[0];
        j_old = RC->Gc->data[ix];
        i7 = (int32_T)emlrtIntegerCheckFastR2012b(j_old, &k_emlrtDCI,
          emlrtRootTLSGlobal);
        emlrtDynamicBoundsCheckFastR2012b(i7, 1, i6, &k_emlrtBCI,
          emlrtRootTLSGlobal);
      }

      ix = r1->size[0];
      r1->size[0] = v2->size[0];
      emxEnsureCapacity((emxArray__common *)r1, ix, (int32_T)sizeof(boolean_T),
                        &e_emlrtRTEI);
      b_i = v2->size[0];
      for (ix = 0; ix < b_i; ix++) {
        r1->data[ix] = (v2->data[ix] == 0);
      }

      ix = RC->Gc->size[0];
      i6 = r1->size[0];
      emlrtSizeEqCheck1DFastR2012b(ix, i6, &e_emlrtECI, emlrtRootTLSGlobal);
      b_i = RC->Gr->size[0];
      for (ix = 0; ix < b_i; ix++) {
        i6 = dVg->size[0];
        j_old = RC->Gr->data[ix];
        i7 = (int32_T)emlrtIntegerCheckFastR2012b(j_old, &l_emlrtDCI,
          emlrtRootTLSGlobal);
        emlrtDynamicBoundsCheckFastR2012b(i7, 1, i6, &l_emlrtBCI,
          emlrtRootTLSGlobal);
      }

      ix = RC->Gr->size[0];
      i6 = v2->size[0];
      emlrtSizeEqCheck1DFastR2012b(ix, i6, &f_emlrtECI, emlrtRootTLSGlobal);
      ix = r4->size[0];
      r4->size[0] = RC->Gc->size[0];
      emxEnsureCapacity((emxArray__common *)r4, ix, (int32_T)sizeof(real_T),
                        &e_emlrtRTEI);
      b_i = RC->Gc->size[0];
      for (ix = 0; ix < b_i; ix++) {
        r4->data[ix] = dVg->data[(int32_T)RC->Gc->data[ix] - 1] * (real_T)
          r1->data[ix];
      }

      ix = r5->size[0];
      r5->size[0] = RC->Gr->size[0];
      emxEnsureCapacity((emxArray__common *)r5, ix, (int32_T)sizeof(real_T),
                        &e_emlrtRTEI);
      b_i = RC->Gr->size[0];
      for (ix = 0; ix < b_i; ix++) {
        r5->data[ix] = dVg->data[(int32_T)RC->Gr->data[ix] - 1] * (real_T)
          v2->data[ix];
      }

      ix = r4->size[0];
      i6 = r5->size[0];
      emlrtSizeEqCheck1DFastR2012b(ix, i6, &e_emlrtECI, emlrtRootTLSGlobal);
      ix = CL->size[0];
      i6 = dPc->size[0];
      emlrtSizeEqCheck1DFastR2012b(ix, i6, &d_emlrtECI, emlrtRootTLSGlobal);
      ix = GL->size[0];
      i6 = dPg->size[0];
      emlrtSizeEqCheck1DFastR2012b(ix, i6, &c_emlrtECI, emlrtRootTLSGlobal);
      emlrtPushRtStackR2012b(&c_emlrtRSI, emlrtRootTLSGlobal);
      ix = b_GL->size[0];
      b_GL->size[0] = GL->size[0];
      emxEnsureCapacity((emxArray__common *)b_GL, ix, (int32_T)sizeof(real_T),
                        &e_emlrtRTEI);
      b_i = GL->size[0];
      for (ix = 0; ix < b_i; ix++) {
        b_GL->data[ix] = GL->data[ix] * dPg->data[ix];
      }

      ix = r10->size[0];
      r10->size[0] = r4->size[0];
      emxEnsureCapacity((emxArray__common *)r10, ix, (int32_T)sizeof(real_T),
                        &e_emlrtRTEI);
      b_i = r4->size[0];
      for (ix = 0; ix < b_i; ix++) {
        r10->data[ix] = r4->data[ix] + r5->data[ix];
      }

      rdivide(b_GL, r10, r4);
      ix = b_CL->size[0];
      b_CL->size[0] = CL->size[0];
      emxEnsureCapacity((emxArray__common *)b_CL, ix, (int32_T)sizeof(real_T),
                        &e_emlrtRTEI);
      b_i = CL->size[0];
      for (ix = 0; ix < b_i; ix++) {
        b_CL->data[ix] = CL->data[ix] * dPc->data[ix];
      }

      ix = r9->size[0];
      r9->size[0] = r2->size[0];
      emxEnsureCapacity((emxArray__common *)r9, ix, (int32_T)sizeof(real_T),
                        &e_emlrtRTEI);
      b_i = r2->size[0];
      for (ix = 0; ix < b_i; ix++) {
        r9->data[ix] = r2->data[ix] + r3->data[ix];
      }

      rdivide(b_CL, r9, r2);
      ix = r8->size[0];
      r8->size[0] = r2->size[0] + r4->size[0];
      emxEnsureCapacity((emxArray__common *)r8, ix, (int32_T)sizeof(real_T),
                        &e_emlrtRTEI);
      b_i = r2->size[0];
      for (ix = 0; ix < b_i; ix++) {
        r8->data[ix] = r2->data[ix];
      }

      b_i = r4->size[0];
      for (ix = 0; ix < b_i; ix++) {
        r8->data[ix + r2->size[0]] = r4->data[ix];
      }

      b_abs(r8, dPg);
      emlrtPushRtStackR2012b(&r_emlrtRSI, emlrtRootTLSGlobal);
      emlrtPushRtStackR2012b(&f_emlrtRSI, emlrtRootTLSGlobal);
      if ((dPg->size[0] == 1) || (dPg->size[0] != 1)) {
        overflow = TRUE;
      } else {
        overflow = FALSE;
      }

      if (overflow) {
      } else {
        emlrtPushRtStackR2012b(&g_emlrtRSI, emlrtRootTLSGlobal);
        e_y = NULL;
        m8 = mxCreateCharArray(2, iv19);
        for (b_i = 0; b_i < 36; b_i++) {
          cv6[b_i] = cv7[b_i];
        }

        emlrtInitCharArrayR2013a(emlrtRootTLSGlobal, 36, m8, cv6);
        emlrtAssign(&e_y, m8);
        error(message(e_y, &emlrtMCI), &b_emlrtMCI);
        emlrtPopRtStackR2012b(&g_emlrtRSI, emlrtRootTLSGlobal);
      }

      if (dPg->size[0] > 0) {
      } else {
        emlrtPushRtStackR2012b(&h_emlrtRSI, emlrtRootTLSGlobal);
        f_y = NULL;
        m8 = mxCreateCharArray(2, iv20);
        for (b_i = 0; b_i < 39; b_i++) {
          cv8[b_i] = cv9[b_i];
        }

        emlrtInitCharArrayR2013a(emlrtRootTLSGlobal, 39, m8, cv8);
        emlrtAssign(&f_y, m8);
        error(message(f_y, &c_emlrtMCI), &d_emlrtMCI);
        emlrtPopRtStackR2012b(&h_emlrtRSI, emlrtRootTLSGlobal);
      }

      emlrtPushRtStackR2012b(&i_emlrtRSI, emlrtRootTLSGlobal);
      b_i = 1;
      mtmp = dPg->data[0];
      if (dPg->size[0] > 1) {
        if (muDoubleScalarIsNaN(dPg->data[0])) {
          emlrtPushRtStackR2012b(&k_emlrtRSI, emlrtRootTLSGlobal);
          if (2 > dPg->size[0]) {
            overflow = FALSE;
          } else {
            overflow = (dPg->size[0] > 2147483646);
          }

          if (overflow) {
            emlrtPushRtStackR2012b(&l_emlrtRSI, emlrtRootTLSGlobal);
            check_forloop_overflow_error();
            emlrtPopRtStackR2012b(&l_emlrtRSI, emlrtRootTLSGlobal);
          }

          emlrtPopRtStackR2012b(&k_emlrtRSI, emlrtRootTLSGlobal);
          ix = 2;
          exitg3 = FALSE;
          while ((exitg3 == FALSE) && (ix <= dPg->size[0])) {
            b_i = ix;
            if (!muDoubleScalarIsNaN(dPg->data[ix - 1])) {
              mtmp = dPg->data[ix - 1];
              exitg3 = TRUE;
            } else {
              ix++;
            }
          }
        }

        if (b_i < dPg->size[0]) {
          emlrtPushRtStackR2012b(&j_emlrtRSI, emlrtRootTLSGlobal);
          if (b_i + 1 > dPg->size[0]) {
            overflow = FALSE;
          } else {
            overflow = (dPg->size[0] > 2147483646);
          }

          if (overflow) {
            emlrtPushRtStackR2012b(&l_emlrtRSI, emlrtRootTLSGlobal);
            check_forloop_overflow_error();
            emlrtPopRtStackR2012b(&l_emlrtRSI, emlrtRootTLSGlobal);
          }

          emlrtPopRtStackR2012b(&j_emlrtRSI, emlrtRootTLSGlobal);
          while (b_i + 1 <= dPg->size[0]) {
            if (dPg->data[b_i] > mtmp) {
              mtmp = dPg->data[b_i];
            }

            b_i++;
          }
        }
      }

      emlrtPopRtStackR2012b(&i_emlrtRSI, emlrtRootTLSGlobal);
      emlrtPopRtStackR2012b(&f_emlrtRSI, emlrtRootTLSGlobal);
      emlrtPopRtStackR2012b(&r_emlrtRSI, emlrtRootTLSGlobal);
      emlrtPopRtStackR2012b(&c_emlrtRSI, emlrtRootTLSGlobal);
      ix = W1C->size[0];
      i6 = dPwc->size[0];
      emlrtSizeEqCheck1DFastR2012b(ix, i6, &b_emlrtECI, emlrtRootTLSGlobal);
      ix = W1G->size[0];
      i6 = dPwg->size[0];
      emlrtSizeEqCheck1DFastR2012b(ix, i6, &emlrtECI, emlrtRootTLSGlobal);
      emlrtPushRtStackR2012b(&d_emlrtRSI, emlrtRootTLSGlobal);
      ix = b_W1G->size[0];
      b_W1G->size[0] = W1G->size[0];
      emxEnsureCapacity((emxArray__common *)b_W1G, ix, (int32_T)sizeof(real_T),
                        &e_emlrtRTEI);
      b_i = W1G->size[0];
      for (ix = 0; ix < b_i; ix++) {
        b_W1G->data[ix] = W1G->data[ix] * dPwg->data[ix];
      }

      b_i = WoG->size[0];
      ix = b_dVg->size[0];
      b_dVg->size[0] = b_i;
      emxEnsureCapacity((emxArray__common *)b_dVg, ix, (int32_T)sizeof(real_T),
                        &e_emlrtRTEI);
      for (ix = 0; ix < b_i; ix++) {
        i6 = dVg->size[0];
        j_old = WoG->data[ix];
        i7 = (int32_T)emlrtIntegerCheckFastR2012b(j_old, &m_emlrtDCI,
          emlrtRootTLSGlobal);
        b_dVg->data[ix] = dVg->data[emlrtDynamicBoundsCheckFastR2012b(i7, 1, i6,
          &m_emlrtBCI, emlrtRootTLSGlobal) - 1];
      }

      rdivide(b_W1G, b_dVg, r2);
      ix = b_W1C->size[0];
      b_W1C->size[0] = W1C->size[0];
      emxEnsureCapacity((emxArray__common *)b_W1C, ix, (int32_T)sizeof(real_T),
                        &e_emlrtRTEI);
      b_i = W1C->size[0];
      for (ix = 0; ix < b_i; ix++) {
        b_W1C->data[ix] = W1C->data[ix] * dPwc->data[ix];
      }

      b_i = WoC->size[0];
      ix = b_dVc->size[0];
      b_dVc->size[0] = b_i;
      emxEnsureCapacity((emxArray__common *)b_dVc, ix, (int32_T)sizeof(real_T),
                        &e_emlrtRTEI);
      for (ix = 0; ix < b_i; ix++) {
        i6 = dVc->size[0];
        j_old = WoC->data[ix];
        i7 = (int32_T)emlrtIntegerCheckFastR2012b(j_old, &n_emlrtDCI,
          emlrtRootTLSGlobal);
        b_dVc->data[ix] = dVc->data[emlrtDynamicBoundsCheckFastR2012b(i7, 1, i6,
          &n_emlrtBCI, emlrtRootTLSGlobal) - 1];
      }

      rdivide(b_W1C, b_dVc, r3);
      ix = r7->size[0];
      r7->size[0] = r3->size[0] + r2->size[0];
      emxEnsureCapacity((emxArray__common *)r7, ix, (int32_T)sizeof(real_T),
                        &e_emlrtRTEI);
      b_i = r3->size[0];
      for (ix = 0; ix < b_i; ix++) {
        r7->data[ix] = r3->data[ix];
      }

      b_i = r2->size[0];
      for (ix = 0; ix < b_i; ix++) {
        r7->data[ix + r3->size[0]] = r2->data[ix];
      }

      b_abs(r7, dPg);
      emlrtPushRtStackR2012b(&r_emlrtRSI, emlrtRootTLSGlobal);
      emlrtPushRtStackR2012b(&f_emlrtRSI, emlrtRootTLSGlobal);
      if ((dPg->size[0] == 1) || (dPg->size[0] != 1)) {
        overflow = TRUE;
      } else {
        overflow = FALSE;
      }

      if (overflow) {
      } else {
        emlrtPushRtStackR2012b(&g_emlrtRSI, emlrtRootTLSGlobal);
        g_y = NULL;
        m8 = mxCreateCharArray(2, iv21);
        for (b_i = 0; b_i < 36; b_i++) {
          cv6[b_i] = cv7[b_i];
        }

        emlrtInitCharArrayR2013a(emlrtRootTLSGlobal, 36, m8, cv6);
        emlrtAssign(&g_y, m8);
        error(message(g_y, &emlrtMCI), &b_emlrtMCI);
        emlrtPopRtStackR2012b(&g_emlrtRSI, emlrtRootTLSGlobal);
      }

      if (dPg->size[0] > 0) {
      } else {
        emlrtPushRtStackR2012b(&h_emlrtRSI, emlrtRootTLSGlobal);
        h_y = NULL;
        m8 = mxCreateCharArray(2, iv22);
        for (b_i = 0; b_i < 39; b_i++) {
          cv8[b_i] = cv9[b_i];
        }

        emlrtInitCharArrayR2013a(emlrtRootTLSGlobal, 39, m8, cv8);
        emlrtAssign(&h_y, m8);
        error(message(h_y, &c_emlrtMCI), &d_emlrtMCI);
        emlrtPopRtStackR2012b(&h_emlrtRSI, emlrtRootTLSGlobal);
      }

      emlrtPushRtStackR2012b(&i_emlrtRSI, emlrtRootTLSGlobal);
      b_i = 1;
      j_old = dPg->data[0];
      if (dPg->size[0] > 1) {
        if (muDoubleScalarIsNaN(dPg->data[0])) {
          emlrtPushRtStackR2012b(&k_emlrtRSI, emlrtRootTLSGlobal);
          if (2 > dPg->size[0]) {
            overflow = FALSE;
          } else {
            overflow = (dPg->size[0] > 2147483646);
          }

          if (overflow) {
            emlrtPushRtStackR2012b(&l_emlrtRSI, emlrtRootTLSGlobal);
            check_forloop_overflow_error();
            emlrtPopRtStackR2012b(&l_emlrtRSI, emlrtRootTLSGlobal);
          }

          emlrtPopRtStackR2012b(&k_emlrtRSI, emlrtRootTLSGlobal);
          ix = 2;
          exitg2 = FALSE;
          while ((exitg2 == FALSE) && (ix <= dPg->size[0])) {
            b_i = ix;
            if (!muDoubleScalarIsNaN(dPg->data[ix - 1])) {
              j_old = dPg->data[ix - 1];
              exitg2 = TRUE;
            } else {
              ix++;
            }
          }
        }

        if (b_i < dPg->size[0]) {
          emlrtPushRtStackR2012b(&j_emlrtRSI, emlrtRootTLSGlobal);
          if (b_i + 1 > dPg->size[0]) {
            overflow = FALSE;
          } else {
            overflow = (dPg->size[0] > 2147483646);
          }

          if (overflow) {
            emlrtPushRtStackR2012b(&l_emlrtRSI, emlrtRootTLSGlobal);
            check_forloop_overflow_error();
            emlrtPopRtStackR2012b(&l_emlrtRSI, emlrtRootTLSGlobal);
          }

          emlrtPopRtStackR2012b(&j_emlrtRSI, emlrtRootTLSGlobal);
          while (b_i + 1 <= dPg->size[0]) {
            if (dPg->data[b_i] > j_old) {
              j_old = dPg->data[b_i];
            }

            b_i++;
          }
        }
      }

      emlrtPopRtStackR2012b(&i_emlrtRSI, emlrtRootTLSGlobal);
      emlrtPopRtStackR2012b(&f_emlrtRSI, emlrtRootTLSGlobal);
      emlrtPopRtStackR2012b(&r_emlrtRSI, emlrtRootTLSGlobal);
      emlrtPopRtStackR2012b(&d_emlrtRSI, emlrtRootTLSGlobal);
      mtmp = dt / (1.0 / mtmp);
      j_old = dt / (1.0 / j_old);
      b_i = 1;
      if (muDoubleScalarIsNaN(mtmp)) {
        ix = 2;
        exitg1 = FALSE;
        while ((exitg1 == FALSE) && (ix < 3)) {
          b_i = 2;
          if (!muDoubleScalarIsNaN(j_old)) {
            mtmp = j_old;
            exitg1 = TRUE;
          } else {
            ix = 3;
          }
        }
      }

      if ((b_i < 2) && (j_old > mtmp)) {
        mtmp = j_old;
      }

      *ndt = muDoubleScalarCeil(mtmp * PR->Fc2);
    }

    emxFree_real_T(&b_Sw2);
    emxFree_boolean_T(&r11);
    emxFree_real_T(&r10);
    emxFree_real_T(&b_GL);
    emxFree_real_T(&r9);
    emxFree_real_T(&b_CL);
    emxFree_real_T(&r8);
    emxFree_real_T(&b_dVg);
    emxFree_real_T(&b_W1G);
    emxFree_real_T(&b_dVc);
    emxFree_real_T(&b_W1C);
    emxFree_real_T(&r7);
    emxFree_int32_T(&r6);
    emxFree_real_T(&r5);
    emxFree_real_T(&r4);
    emxFree_real_T(&r3);
    emxFree_real_T(&r2);
    emxFree_boolean_T(&r1);
  } else {
    *ndt = 1.0;
  }

  emxFree_boolean_T(&v2);
  emxFree_boolean_T(&v1);
  emxFree_real_T(&dPwg);
  emxFree_real_T(&dPwc);
  emxFree_real_T(&dPg);
  emxFree_real_T(&dPc);
  j_old = *j_ndt;
  *j_ndt += 1.0 / *ndt;
  fl1 = (*j_ndt < 1.0);
  if (fl1 == 0) {
    *ndt = 1.0 / (1.0 - j_old);
  }

  emlrtHeapReferenceStackLeaveFcnR2012b(emlrtRootTLSGlobal);
  return fl1;
}

/* End of code generation (vibor_t.c) */
