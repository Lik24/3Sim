/*
 * vibor_t_api.c
 *
 * Code generation for function 'vibor_t_api'
 *
 * C source code generated on: Tue Aug 12 14:07:28 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "vibor_t.h"
#include "vibor_t_api.h"
#include "vibor_t_emxutil.h"

/* Variable Definitions */
static emlrtRTEInfo d_emlrtRTEI = { 1, 1, "vibor_t_api", "" };

/* Function Declarations */
static void ab_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier
  *msgId, emxArray_real_T *ret);
static real_T b_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier
  *parentId);
static const mxArray *b_emlrt_marshallOut(boolean_T u);
static void bb_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier
  *msgId, emxArray_real_T *ret);
static void c_emlrt_marshallIn(const mxArray *Pc, const char_T *identifier,
  emxArray_real_T *y);
static void cb_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier
  *msgId, emxArray_real_T *ret);
static void d_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier
  *parentId, emxArray_real_T *y);
static void e_emlrt_marshallIn(const mxArray *PR, const char_T *identifier,
  b_struct_T *y);
static real_T emlrt_marshallIn(const mxArray *ndt, const char_T *identifier);
static const mxArray *emlrt_marshallOut(real_T u);
static void f_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier
  *parentId, b_struct_T *y);
static void g_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier
  *parentId, real_T y[3]);
static void h_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier
  *parentId, real_T y[4]);
static void i_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier
  *parentId, real_T y[5]);
static void info_helper(ResolvedFunctionInfo info[35]);
static void j_emlrt_marshallIn(const mxArray *RC, const char_T *identifier,
  struct_T *y);
static void k_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier
  *parentId, struct_T *y);
static void l_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier
  *parentId, emxArray_real_T *y);
static void m_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier
  *parentId, emxArray_real_T *y);
static void n_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier
  *parentId, emxArray_real_T *y);
static void o_emlrt_marshallIn(const mxArray *GL, const char_T *identifier,
  emxArray_real_T *y);
static void p_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier
  *parentId, emxArray_real_T *y);
static void q_emlrt_marshallIn(const mxArray *WoC, const char_T *identifier,
  emxArray_real_T *y);
static void r_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier
  *parentId, emxArray_real_T *y);
static real_T s_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier
  *msgId);
static void t_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier
  *msgId, emxArray_real_T *ret);
static void u_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier
  *msgId, real_T ret[3]);
static void v_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier
  *msgId, real_T ret[4]);
static void w_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier
  *msgId, real_T ret[5]);
static void x_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier
  *msgId, emxArray_real_T *ret);
static void y_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier
  *msgId, emxArray_real_T *ret);

/* Function Definitions */
static void ab_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier
  *msgId, emxArray_real_T *ret)
{
  int32_T iv13[2];
  boolean_T bv4[2];
  int32_T i5;
  int32_T iv14[2];
  for (i5 = 0; i5 < 2; i5++) {
    iv13[i5] = (i5 << 1) - 1;
    bv4[i5] = TRUE;
  }

  emlrtCheckVsBuiltInR2012b(emlrtRootTLSGlobal, msgId, src, "double", FALSE, 2U,
    iv13, bv4, iv14);
  i5 = ret->size[0] * ret->size[1];
  ret->size[0] = iv14[0];
  ret->size[1] = iv14[1];
  emxEnsureCapacity((emxArray__common *)ret, i5, (int32_T)sizeof(real_T),
                    (emlrtRTEInfo *)NULL);
  emlrtImportArrayR2011b(src, ret->data, 8, FALSE);
  emlrtDestroyArray(&src);
}

static real_T b_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier
  *parentId)
{
  real_T y;
  y = s_emlrt_marshallIn(emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static const mxArray *b_emlrt_marshallOut(boolean_T u)
{
  const mxArray *y;
  const mxArray *m5;
  y = NULL;
  m5 = mxCreateLogicalScalar(u);
  emlrtAssign(&y, m5);
  return y;
}

static void bb_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier
  *msgId, emxArray_real_T *ret)
{
  int32_T iv15[2];
  boolean_T bv5[2];
  int32_T i6;
  int32_T iv16[2];
  for (i6 = 0; i6 < 2; i6++) {
    iv15[i6] = (i6 << 1) - 1;
    bv5[i6] = TRUE;
  }

  emlrtCheckVsBuiltInR2012b(emlrtRootTLSGlobal, msgId, src, "double", FALSE, 2U,
    iv15, bv5, iv16);
  ret->size[0] = iv16[0];
  ret->size[1] = iv16[1];
  ret->allocatedSize = ret->size[0] * ret->size[1];
  ret->data = (real_T *)mxGetData(src);
  ret->canFreeData = FALSE;
  emlrtDestroyArray(&src);
}

static void c_emlrt_marshallIn(const mxArray *Pc, const char_T *identifier,
  emxArray_real_T *y)
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = identifier;
  thisId.fParent = NULL;
  d_emlrt_marshallIn(emlrtAlias(Pc), &thisId, y);
  emlrtDestroyArray(&Pc);
}

static void cb_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier
  *msgId, emxArray_real_T *ret)
{
  int32_T iv17[2];
  boolean_T bv6[2];
  int32_T i7;
  static const boolean_T bv7[2] = { TRUE, FALSE };

  int32_T iv18[2];
  for (i7 = 0; i7 < 2; i7++) {
    iv17[i7] = (i7 << 2) - 1;
    bv6[i7] = bv7[i7];
  }

  emlrtCheckVsBuiltInR2012b(emlrtRootTLSGlobal, msgId, src, "double", FALSE, 2U,
    iv17, bv6, iv18);
  ret->size[0] = iv18[0];
  ret->size[1] = iv18[1];
  ret->allocatedSize = ret->size[0] * ret->size[1];
  ret->data = (real_T *)mxGetData(src);
  ret->canFreeData = FALSE;
  emlrtDestroyArray(&src);
}

static void d_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier
  *parentId, emxArray_real_T *y)
{
  t_emlrt_marshallIn(emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static void e_emlrt_marshallIn(const mxArray *PR, const char_T *identifier,
  b_struct_T *y)
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = identifier;
  thisId.fParent = NULL;
  f_emlrt_marshallIn(emlrtAlias(PR), &thisId, y);
  emlrtDestroyArray(&PR);
}

static real_T emlrt_marshallIn(const mxArray *ndt, const char_T *identifier)
{
  real_T y;
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = identifier;
  thisId.fParent = NULL;
  y = b_emlrt_marshallIn(emlrtAlias(ndt), &thisId);
  emlrtDestroyArray(&ndt);
  return y;
}

static const mxArray *emlrt_marshallOut(real_T u)
{
  const mxArray *y;
  const mxArray *m4;
  y = NULL;
  m4 = mxCreateDoubleScalar(u);
  emlrtAssign(&y, m4);
  return y;
}

static void f_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier
  *parentId, b_struct_T *y)
{
  emlrtMsgIdentifier thisId;
  static const char * fieldNames[22] = { "Nl", "as", "aw", "ts", "tw", "mu",
    "Ns", "Ta", "dt", "ndt", "Ro", "bet", "lam", "Cp", "dh", "Kc", "zc", "Bo",
    "Fc", "Fc2", "Sc2", "dl" };

  thisId.fParent = parentId;
  emlrtCheckStructR2012b(emlrtRootTLSGlobal, parentId, u, 22, fieldNames, 0U, 0);
  thisId.fIdentifier = "Nl";
  y->Nl = b_emlrt_marshallIn(emlrtAlias(emlrtGetFieldR2013a(emlrtRootTLSGlobal,
    u, 0, "Nl")), &thisId);
  thisId.fIdentifier = "as";
  g_emlrt_marshallIn(emlrtAlias(emlrtGetFieldR2013a(emlrtRootTLSGlobal, u, 0,
    "as")), &thisId, y->as);
  thisId.fIdentifier = "aw";
  g_emlrt_marshallIn(emlrtAlias(emlrtGetFieldR2013a(emlrtRootTLSGlobal, u, 0,
    "aw")), &thisId, y->aw);
  thisId.fIdentifier = "ts";
  g_emlrt_marshallIn(emlrtAlias(emlrtGetFieldR2013a(emlrtRootTLSGlobal, u, 0,
    "ts")), &thisId, y->ts);
  thisId.fIdentifier = "tw";
  g_emlrt_marshallIn(emlrtAlias(emlrtGetFieldR2013a(emlrtRootTLSGlobal, u, 0,
    "tw")), &thisId, y->tw);
  thisId.fIdentifier = "mu";
  h_emlrt_marshallIn(emlrtAlias(emlrtGetFieldR2013a(emlrtRootTLSGlobal, u, 0,
    "mu")), &thisId, y->mu);
  thisId.fIdentifier = "Ns";
  y->Ns = b_emlrt_marshallIn(emlrtAlias(emlrtGetFieldR2013a(emlrtRootTLSGlobal,
    u, 0, "Ns")), &thisId);
  thisId.fIdentifier = "Ta";
  y->Ta = b_emlrt_marshallIn(emlrtAlias(emlrtGetFieldR2013a(emlrtRootTLSGlobal,
    u, 0, "Ta")), &thisId);
  thisId.fIdentifier = "dt";
  y->dt = b_emlrt_marshallIn(emlrtAlias(emlrtGetFieldR2013a(emlrtRootTLSGlobal,
    u, 0, "dt")), &thisId);
  thisId.fIdentifier = "ndt";
  y->ndt = b_emlrt_marshallIn(emlrtAlias(emlrtGetFieldR2013a(emlrtRootTLSGlobal,
    u, 0, "ndt")), &thisId);
  thisId.fIdentifier = "Ro";
  i_emlrt_marshallIn(emlrtAlias(emlrtGetFieldR2013a(emlrtRootTLSGlobal, u, 0,
    "Ro")), &thisId, y->Ro);
  thisId.fIdentifier = "bet";
  y->bet = b_emlrt_marshallIn(emlrtAlias(emlrtGetFieldR2013a(emlrtRootTLSGlobal,
    u, 0, "bet")), &thisId);
  thisId.fIdentifier = "lam";
  i_emlrt_marshallIn(emlrtAlias(emlrtGetFieldR2013a(emlrtRootTLSGlobal, u, 0,
    "lam")), &thisId, y->lam);
  thisId.fIdentifier = "Cp";
  i_emlrt_marshallIn(emlrtAlias(emlrtGetFieldR2013a(emlrtRootTLSGlobal, u, 0,
    "Cp")), &thisId, y->Cp);
  thisId.fIdentifier = "dh";
  y->dh = b_emlrt_marshallIn(emlrtAlias(emlrtGetFieldR2013a(emlrtRootTLSGlobal,
    u, 0, "dh")), &thisId);
  thisId.fIdentifier = "Kc";
  y->Kc = b_emlrt_marshallIn(emlrtAlias(emlrtGetFieldR2013a(emlrtRootTLSGlobal,
    u, 0, "Kc")), &thisId);
  thisId.fIdentifier = "zc";
  i_emlrt_marshallIn(emlrtAlias(emlrtGetFieldR2013a(emlrtRootTLSGlobal, u, 0,
    "zc")), &thisId, y->zc);
  thisId.fIdentifier = "Bo";
  i_emlrt_marshallIn(emlrtAlias(emlrtGetFieldR2013a(emlrtRootTLSGlobal, u, 0,
    "Bo")), &thisId, y->Bo);
  thisId.fIdentifier = "Fc";
  y->Fc = b_emlrt_marshallIn(emlrtAlias(emlrtGetFieldR2013a(emlrtRootTLSGlobal,
    u, 0, "Fc")), &thisId);
  thisId.fIdentifier = "Fc2";
  y->Fc2 = b_emlrt_marshallIn(emlrtAlias(emlrtGetFieldR2013a(emlrtRootTLSGlobal,
    u, 0, "Fc2")), &thisId);
  thisId.fIdentifier = "Sc2";
  y->Sc2 = b_emlrt_marshallIn(emlrtAlias(emlrtGetFieldR2013a(emlrtRootTLSGlobal,
    u, 0, "Sc2")), &thisId);
  thisId.fIdentifier = "dl";
  y->dl = b_emlrt_marshallIn(emlrtAlias(emlrtGetFieldR2013a(emlrtRootTLSGlobal,
    u, 0, "dl")), &thisId);
  emlrtDestroyArray(&u);
}

static void g_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier
  *parentId, real_T y[3])
{
  u_emlrt_marshallIn(emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static void h_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier
  *parentId, real_T y[4])
{
  v_emlrt_marshallIn(emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static void i_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier
  *parentId, real_T y[5])
{
  w_emlrt_marshallIn(emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static void info_helper(ResolvedFunctionInfo info[35])
{
  info[0].context =
    "[E]C:/Users/LIK/Documents/MATLAB/Tube_dev/2be_0.661/Sim_Lib/vibor_t.m";
  info[0].name = "min";
  info[0].dominantType = "double";
  info[0].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/min.m";
  info[0].fileTimeLo = 1311240918U;
  info[0].fileTimeHi = 0U;
  info[0].mFileTimeLo = 0U;
  info[0].mFileTimeHi = 0U;
  info[1].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/min.m";
  info[1].name = "eml_min_or_max";
  info[1].dominantType = "char";
  info[1].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m";
  info[1].fileTimeLo = 1334057090U;
  info[1].fileTimeHi = 0U;
  info[1].mFileTimeLo = 0U;
  info[1].mFileTimeHi = 0U;
  info[2].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_extremum";
  info[2].name = "eml_const_nonsingleton_dim";
  info[2].dominantType = "double";
  info[2].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_const_nonsingleton_dim.m";
  info[2].fileTimeLo = 1286804296U;
  info[2].fileTimeHi = 0U;
  info[2].mFileTimeLo = 0U;
  info[2].mFileTimeHi = 0U;
  info[3].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_extremum";
  info[3].name = "eml_scalar_eg";
  info[3].dominantType = "double";
  info[3].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m";
  info[3].fileTimeLo = 1286804396U;
  info[3].fileTimeHi = 0U;
  info[3].mFileTimeLo = 0U;
  info[3].mFileTimeHi = 0U;
  info[4].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_extremum";
  info[4].name = "eml_index_class";
  info[4].dominantType = "";
  info[4].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m";
  info[4].fileTimeLo = 1323156178U;
  info[4].fileTimeHi = 0U;
  info[4].mFileTimeLo = 0U;
  info[4].mFileTimeHi = 0U;
  info[5].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_extremum_sub";
  info[5].name = "eml_index_class";
  info[5].dominantType = "";
  info[5].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m";
  info[5].fileTimeLo = 1323156178U;
  info[5].fileTimeHi = 0U;
  info[5].mFileTimeLo = 0U;
  info[5].mFileTimeHi = 0U;
  info[6].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_extremum_sub";
  info[6].name = "isnan";
  info[6].dominantType = "double";
  info[6].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isnan.m";
  info[6].fileTimeLo = 1286804360U;
  info[6].fileTimeHi = 0U;
  info[6].mFileTimeLo = 0U;
  info[6].mFileTimeHi = 0U;
  info[7].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_extremum_sub";
  info[7].name = "eml_index_plus";
  info[7].dominantType = "coder.internal.indexInt";
  info[7].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m";
  info[7].fileTimeLo = 1286804378U;
  info[7].fileTimeHi = 0U;
  info[7].mFileTimeLo = 0U;
  info[7].mFileTimeHi = 0U;
  info[8].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m";
  info[8].name = "eml_index_class";
  info[8].dominantType = "";
  info[8].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m";
  info[8].fileTimeLo = 1323156178U;
  info[8].fileTimeHi = 0U;
  info[8].mFileTimeLo = 0U;
  info[8].mFileTimeHi = 0U;
  info[9].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_extremum_sub";
  info[9].name = "eml_int_forloop_overflow_check";
  info[9].dominantType = "";
  info[9].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m";
  info[9].fileTimeLo = 1346495940U;
  info[9].fileTimeHi = 0U;
  info[9].mFileTimeLo = 0U;
  info[9].mFileTimeHi = 0U;
  info[10].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m!eml_int_forloop_overflow_check_helper";
  info[10].name = "intmax";
  info[10].dominantType = "char";
  info[10].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m";
  info[10].fileTimeLo = 1311240916U;
  info[10].fileTimeHi = 0U;
  info[10].mFileTimeLo = 0U;
  info[10].mFileTimeHi = 0U;
  info[11].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_extremum_sub";
  info[11].name = "eml_relop";
  info[11].dominantType = "function_handle";
  info[11].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_relop.m";
  info[11].fileTimeLo = 1342436782U;
  info[11].fileTimeHi = 0U;
  info[11].mFileTimeLo = 0U;
  info[11].mFileTimeHi = 0U;
  info[12].context =
    "[E]C:/Users/LIK/Documents/MATLAB/Tube_dev/2be_0.661/Sim_Lib/vibor_t.m";
  info[12].name = "abs";
  info[12].dominantType = "double";
  info[12].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m";
  info[12].fileTimeLo = 1343815966U;
  info[12].fileTimeHi = 0U;
  info[12].mFileTimeLo = 0U;
  info[12].mFileTimeHi = 0U;
  info[13].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m";
  info[13].name = "eml_scalar_abs";
  info[13].dominantType = "double";
  info[13].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_abs.m";
  info[13].fileTimeLo = 1286804312U;
  info[13].fileTimeHi = 0U;
  info[13].mFileTimeLo = 0U;
  info[13].mFileTimeHi = 0U;
  info[14].context =
    "[E]C:/Users/LIK/Documents/MATLAB/Tube_dev/2be_0.661/Sim_Lib/vibor_t.m";
  info[14].name = "isnan";
  info[14].dominantType = "double";
  info[14].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isnan.m";
  info[14].fileTimeLo = 1286804360U;
  info[14].fileTimeHi = 0U;
  info[14].mFileTimeLo = 0U;
  info[14].mFileTimeHi = 0U;
  info[15].context =
    "[E]C:/Users/LIK/Documents/MATLAB/Tube_dev/2be_0.661/Sim_Lib/vibor_t.m";
  info[15].name = "eml_li_find";
  info[15].dominantType = "";
  info[15].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_li_find.m";
  info[15].fileTimeLo = 1286804386U;
  info[15].fileTimeHi = 0U;
  info[15].mFileTimeLo = 0U;
  info[15].mFileTimeHi = 0U;
  info[16].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_li_find.m";
  info[16].name = "eml_index_class";
  info[16].dominantType = "";
  info[16].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m";
  info[16].fileTimeLo = 1323156178U;
  info[16].fileTimeHi = 0U;
  info[16].mFileTimeLo = 0U;
  info[16].mFileTimeHi = 0U;
  info[17].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_li_find.m!compute_nones";
  info[17].name = "eml_index_class";
  info[17].dominantType = "";
  info[17].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m";
  info[17].fileTimeLo = 1323156178U;
  info[17].fileTimeHi = 0U;
  info[17].mFileTimeLo = 0U;
  info[17].mFileTimeHi = 0U;
  info[18].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_li_find.m!compute_nones";
  info[18].name = "eml_int_forloop_overflow_check";
  info[18].dominantType = "";
  info[18].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m";
  info[18].fileTimeLo = 1346495940U;
  info[18].fileTimeHi = 0U;
  info[18].mFileTimeLo = 0U;
  info[18].mFileTimeHi = 0U;
  info[19].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_li_find.m!compute_nones";
  info[19].name = "eml_index_plus";
  info[19].dominantType = "double";
  info[19].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m";
  info[19].fileTimeLo = 1286804378U;
  info[19].fileTimeHi = 0U;
  info[19].mFileTimeLo = 0U;
  info[19].mFileTimeHi = 0U;
  info[20].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_li_find.m";
  info[20].name = "eml_int_forloop_overflow_check";
  info[20].dominantType = "";
  info[20].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m";
  info[20].fileTimeLo = 1346495940U;
  info[20].fileTimeHi = 0U;
  info[20].mFileTimeLo = 0U;
  info[20].mFileTimeHi = 0U;
  info[21].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_li_find.m";
  info[21].name = "eml_index_plus";
  info[21].dominantType = "double";
  info[21].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m";
  info[21].fileTimeLo = 1286804378U;
  info[21].fileTimeHi = 0U;
  info[21].mFileTimeLo = 0U;
  info[21].mFileTimeHi = 0U;
  info[22].context =
    "[E]C:/Users/LIK/Documents/MATLAB/Tube_dev/2be_0.661/Sim_Lib/vibor_t.m";
  info[22].name = "max";
  info[22].dominantType = "double";
  info[22].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/max.m";
  info[22].fileTimeLo = 1311240916U;
  info[22].fileTimeHi = 0U;
  info[22].mFileTimeLo = 0U;
  info[22].mFileTimeHi = 0U;
  info[23].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/max.m";
  info[23].name = "eml_min_or_max";
  info[23].dominantType = "char";
  info[23].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m";
  info[23].fileTimeLo = 1334057090U;
  info[23].fileTimeHi = 0U;
  info[23].mFileTimeLo = 0U;
  info[23].mFileTimeHi = 0U;
  info[24].context =
    "[E]C:/Users/LIK/Documents/MATLAB/Tube_dev/2be_0.661/Sim_Lib/vibor_t.m";
  info[24].name = "mrdivide";
  info[24].dominantType = "double";
  info[24].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p";
  info[24].fileTimeLo = 1357937148U;
  info[24].fileTimeHi = 0U;
  info[24].mFileTimeLo = 1319715566U;
  info[24].mFileTimeHi = 0U;
  info[25].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p";
  info[25].name = "rdivide";
  info[25].dominantType = "double";
  info[25].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m";
  info[25].fileTimeLo = 1346495988U;
  info[25].fileTimeHi = 0U;
  info[25].mFileTimeLo = 0U;
  info[25].mFileTimeHi = 0U;
  info[26].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m";
  info[26].name = "eml_scalexp_compatible";
  info[26].dominantType = "double";
  info[26].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_compatible.m";
  info[26].fileTimeLo = 1286804396U;
  info[26].fileTimeHi = 0U;
  info[26].mFileTimeLo = 0U;
  info[26].mFileTimeHi = 0U;
  info[27].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m";
  info[27].name = "eml_div";
  info[27].dominantType = "double";
  info[27].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m";
  info[27].fileTimeLo = 1313333410U;
  info[27].fileTimeHi = 0U;
  info[27].mFileTimeLo = 0U;
  info[27].mFileTimeHi = 0U;
  info[28].context =
    "[E]C:/Users/LIK/Documents/MATLAB/Tube_dev/2be_0.661/Sim_Lib/vibor_t.m";
  info[28].name = "mtimes";
  info[28].dominantType = "double";
  info[28].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m";
  info[28].fileTimeLo = 1289505292U;
  info[28].fileTimeHi = 0U;
  info[28].mFileTimeLo = 0U;
  info[28].mFileTimeHi = 0U;
  info[29].context =
    "[E]C:/Users/LIK/Documents/MATLAB/Tube_dev/2be_0.661/Sim_Lib/vibor_t.m";
  info[29].name = "ceil";
  info[29].dominantType = "double";
  info[29].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/ceil.m";
  info[29].fileTimeLo = 1343815972U;
  info[29].fileTimeHi = 0U;
  info[29].mFileTimeLo = 0U;
  info[29].mFileTimeHi = 0U;
  info[30].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/ceil.m";
  info[30].name = "eml_scalar_ceil";
  info[30].dominantType = "double";
  info[30].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_ceil.m";
  info[30].fileTimeLo = 1286804320U;
  info[30].fileTimeHi = 0U;
  info[30].mFileTimeLo = 0U;
  info[30].mFileTimeHi = 0U;
  info[31].context =
    "[E]C:/Users/LIK/Documents/MATLAB/Tube_dev/2be_0.661/Sim_Lib/vibor_t.m";
  info[31].name = "rdivide";
  info[31].dominantType = "double";
  info[31].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m";
  info[31].fileTimeLo = 1346495988U;
  info[31].fileTimeHi = 0U;
  info[31].mFileTimeLo = 0U;
  info[31].mFileTimeHi = 0U;
  info[32].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_compatible.m";
  info[32].name = "isequal";
  info[32].dominantType = "double";
  info[32].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isequal.m";
  info[32].fileTimeLo = 1286804358U;
  info[32].fileTimeHi = 0U;
  info[32].mFileTimeLo = 0U;
  info[32].mFileTimeHi = 0U;
  info[33].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isequal.m";
  info[33].name = "eml_isequal_core";
  info[33].dominantType = "double";
  info[33].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_isequal_core.m";
  info[33].fileTimeLo = 1286804386U;
  info[33].fileTimeHi = 0U;
  info[33].mFileTimeLo = 0U;
  info[33].mFileTimeHi = 0U;
  info[34].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_isequal_core.m!isequal_scalar";
  info[34].name = "isnan";
  info[34].dominantType = "double";
  info[34].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isnan.m";
  info[34].fileTimeLo = 1286804360U;
  info[34].fileTimeHi = 0U;
  info[34].mFileTimeLo = 0U;
  info[34].mFileTimeHi = 0U;
}

static void j_emlrt_marshallIn(const mxArray *RC, const char_T *identifier,
  struct_T *y)
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = identifier;
  thisId.fParent = NULL;
  k_emlrt_marshallIn(emlrtAlias(RC), &thisId, y);
  emlrtDestroyArray(&RC);
}

static void k_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier
  *parentId, struct_T *y)
{
  emlrtMsgIdentifier thisId;
  static const char * fieldNames[13] = { "rc", "Arc", "Cr", "Cc", "Gr", "Gc",
    "ACr", "ACc", "AGr", "AGc", "na", "nc", "ng" };

  thisId.fParent = parentId;
  emlrtCheckStructR2012b(emlrtRootTLSGlobal, parentId, u, 13, fieldNames, 0U, 0);
  thisId.fIdentifier = "rc";
  l_emlrt_marshallIn(emlrtAlias(emlrtGetFieldR2013a(emlrtRootTLSGlobal, u, 0,
    "rc")), &thisId, y->rc);
  thisId.fIdentifier = "Arc";
  l_emlrt_marshallIn(emlrtAlias(emlrtGetFieldR2013a(emlrtRootTLSGlobal, u, 0,
    "Arc")), &thisId, y->Arc);
  thisId.fIdentifier = "Cr";
  m_emlrt_marshallIn(emlrtAlias(emlrtGetFieldR2013a(emlrtRootTLSGlobal, u, 0,
    "Cr")), &thisId, y->Cr);
  thisId.fIdentifier = "Cc";
  m_emlrt_marshallIn(emlrtAlias(emlrtGetFieldR2013a(emlrtRootTLSGlobal, u, 0,
    "Cc")), &thisId, y->Cc);
  thisId.fIdentifier = "Gr";
  n_emlrt_marshallIn(emlrtAlias(emlrtGetFieldR2013a(emlrtRootTLSGlobal, u, 0,
    "Gr")), &thisId, y->Gr);
  thisId.fIdentifier = "Gc";
  n_emlrt_marshallIn(emlrtAlias(emlrtGetFieldR2013a(emlrtRootTLSGlobal, u, 0,
    "Gc")), &thisId, y->Gc);
  thisId.fIdentifier = "ACr";
  m_emlrt_marshallIn(emlrtAlias(emlrtGetFieldR2013a(emlrtRootTLSGlobal, u, 0,
    "ACr")), &thisId, y->ACr);
  thisId.fIdentifier = "ACc";
  m_emlrt_marshallIn(emlrtAlias(emlrtGetFieldR2013a(emlrtRootTLSGlobal, u, 0,
    "ACc")), &thisId, y->ACc);
  thisId.fIdentifier = "AGr";
  m_emlrt_marshallIn(emlrtAlias(emlrtGetFieldR2013a(emlrtRootTLSGlobal, u, 0,
    "AGr")), &thisId, y->AGr);
  thisId.fIdentifier = "AGc";
  m_emlrt_marshallIn(emlrtAlias(emlrtGetFieldR2013a(emlrtRootTLSGlobal, u, 0,
    "AGc")), &thisId, y->AGc);
  thisId.fIdentifier = "na";
  y->na = b_emlrt_marshallIn(emlrtAlias(emlrtGetFieldR2013a(emlrtRootTLSGlobal,
    u, 0, "na")), &thisId);
  thisId.fIdentifier = "nc";
  y->nc = b_emlrt_marshallIn(emlrtAlias(emlrtGetFieldR2013a(emlrtRootTLSGlobal,
    u, 0, "nc")), &thisId);
  thisId.fIdentifier = "ng";
  y->ng = b_emlrt_marshallIn(emlrtAlias(emlrtGetFieldR2013a(emlrtRootTLSGlobal,
    u, 0, "ng")), &thisId);
  emlrtDestroyArray(&u);
}

static void l_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier
  *parentId, emxArray_real_T *y)
{
  x_emlrt_marshallIn(emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static void m_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier
  *parentId, emxArray_real_T *y)
{
  y_emlrt_marshallIn(emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static void n_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier
  *parentId, emxArray_real_T *y)
{
  ab_emlrt_marshallIn(emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static void o_emlrt_marshallIn(const mxArray *GL, const char_T *identifier,
  emxArray_real_T *y)
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = identifier;
  thisId.fParent = NULL;
  p_emlrt_marshallIn(emlrtAlias(GL), &thisId, y);
  emlrtDestroyArray(&GL);
}

static void p_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier
  *parentId, emxArray_real_T *y)
{
  bb_emlrt_marshallIn(emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static void q_emlrt_marshallIn(const mxArray *WoC, const char_T *identifier,
  emxArray_real_T *y)
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = identifier;
  thisId.fParent = NULL;
  r_emlrt_marshallIn(emlrtAlias(WoC), &thisId, y);
  emlrtDestroyArray(&WoC);
}

static void r_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier
  *parentId, emxArray_real_T *y)
{
  cb_emlrt_marshallIn(emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static real_T s_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier
  *msgId)
{
  real_T ret;
  emlrtCheckBuiltInR2012b(emlrtRootTLSGlobal, msgId, src, "double", FALSE, 0U, 0);
  ret = *(real_T *)mxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

static void t_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier
  *msgId, emxArray_real_T *ret)
{
  int32_T iv4[1];
  boolean_T bv0[1];
  int32_T iv5[1];
  iv4[0] = -1;
  bv0[0] = TRUE;
  emlrtCheckVsBuiltInR2012b(emlrtRootTLSGlobal, msgId, src, "double", FALSE, 1U,
    iv4, bv0, iv5);
  ret->size[0] = iv5[0];
  ret->allocatedSize = ret->size[0];
  ret->data = (real_T *)mxGetData(src);
  ret->canFreeData = FALSE;
  emlrtDestroyArray(&src);
}

static void u_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier
  *msgId, real_T ret[3])
{
  int32_T iv6[2];
  int32_T i0;
  for (i0 = 0; i0 < 2; i0++) {
    iv6[i0] = 1 + (i0 << 1);
  }

  emlrtCheckBuiltInR2012b(emlrtRootTLSGlobal, msgId, src, "double", FALSE, 2U,
    iv6);
  for (i0 = 0; i0 < 3; i0++) {
    ret[i0] = (*(real_T (*)[3])mxGetData(src))[i0];
  }

  emlrtDestroyArray(&src);
}

static void v_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier
  *msgId, real_T ret[4])
{
  int32_T iv7[2];
  int32_T i1;
  for (i1 = 0; i1 < 2; i1++) {
    iv7[i1] = 1 + 3 * i1;
  }

  emlrtCheckBuiltInR2012b(emlrtRootTLSGlobal, msgId, src, "double", FALSE, 2U,
    iv7);
  for (i1 = 0; i1 < 4; i1++) {
    ret[i1] = (*(real_T (*)[4])mxGetData(src))[i1];
  }

  emlrtDestroyArray(&src);
}

static void w_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier
  *msgId, real_T ret[5])
{
  int32_T iv8[2];
  int32_T i2;
  for (i2 = 0; i2 < 2; i2++) {
    iv8[i2] = 1 + (i2 << 2);
  }

  emlrtCheckBuiltInR2012b(emlrtRootTLSGlobal, msgId, src, "double", FALSE, 2U,
    iv8);
  for (i2 = 0; i2 < 5; i2++) {
    ret[i2] = (*(real_T (*)[5])mxGetData(src))[i2];
  }

  emlrtDestroyArray(&src);
}

static void x_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier
  *msgId, emxArray_real_T *ret)
{
  int32_T iv9[2];
  boolean_T bv1[2];
  int32_T i3;
  static const boolean_T bv2[2] = { TRUE, FALSE };

  int32_T iv10[2];
  for (i3 = 0; i3 < 2; i3++) {
    iv9[i3] = 3 * i3 - 1;
    bv1[i3] = bv2[i3];
  }

  emlrtCheckVsBuiltInR2012b(emlrtRootTLSGlobal, msgId, src, "double", FALSE, 2U,
    iv9, bv1, iv10);
  i3 = ret->size[0] * ret->size[1];
  ret->size[0] = iv10[0];
  ret->size[1] = iv10[1];
  emxEnsureCapacity((emxArray__common *)ret, i3, (int32_T)sizeof(real_T),
                    (emlrtRTEInfo *)NULL);
  emlrtImportArrayR2011b(src, ret->data, 8, FALSE);
  emlrtDestroyArray(&src);
}

static void y_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier
  *msgId, emxArray_real_T *ret)
{
  int32_T iv11[1];
  boolean_T bv3[1];
  int32_T iv12[1];
  int32_T i4;
  iv11[0] = -1;
  bv3[0] = TRUE;
  emlrtCheckVsBuiltInR2012b(emlrtRootTLSGlobal, msgId, src, "double", FALSE, 1U,
    iv11, bv3, iv12);
  i4 = ret->size[0];
  ret->size[0] = iv12[0];
  emxEnsureCapacity((emxArray__common *)ret, i4, (int32_T)sizeof(real_T),
                    (emlrtRTEInfo *)NULL);
  emlrtImportArrayR2011b(src, ret->data, 8, FALSE);
  emlrtDestroyArray(&src);
}

const mxArray *emlrtMexFcnResolvedFunctionsInfo(void)
{
  const mxArray *nameCaptureInfo;
  ResolvedFunctionInfo info[35];
  ResolvedFunctionInfo u[35];
  int32_T i;
  const mxArray *y;
  int32_T iv3[1];
  ResolvedFunctionInfo *r0;
  const char * b_u;
  const mxArray *b_y;
  const mxArray *m3;
  const mxArray *c_y;
  const mxArray *d_y;
  const mxArray *e_y;
  uint32_T c_u;
  const mxArray *f_y;
  const mxArray *g_y;
  const mxArray *h_y;
  const mxArray *i_y;
  nameCaptureInfo = NULL;
  info_helper(info);
  for (i = 0; i < 35; i++) {
    u[i] = info[i];
  }

  y = NULL;
  iv3[0] = 35;
  emlrtAssign(&y, mxCreateStructArray(1, iv3, 0, NULL));
  for (i = 0; i < 35; i++) {
    r0 = &u[i];
    b_u = r0->context;
    b_y = NULL;
    m3 = mxCreateString(b_u);
    emlrtAssign(&b_y, m3);
    emlrtAddField(y, b_y, "context", i);
    b_u = r0->name;
    c_y = NULL;
    m3 = mxCreateString(b_u);
    emlrtAssign(&c_y, m3);
    emlrtAddField(y, c_y, "name", i);
    b_u = r0->dominantType;
    d_y = NULL;
    m3 = mxCreateString(b_u);
    emlrtAssign(&d_y, m3);
    emlrtAddField(y, d_y, "dominantType", i);
    b_u = r0->resolved;
    e_y = NULL;
    m3 = mxCreateString(b_u);
    emlrtAssign(&e_y, m3);
    emlrtAddField(y, e_y, "resolved", i);
    c_u = r0->fileTimeLo;
    f_y = NULL;
    m3 = mxCreateNumericMatrix(1, 1, mxUINT32_CLASS, mxREAL);
    *(uint32_T *)mxGetData(m3) = c_u;
    emlrtAssign(&f_y, m3);
    emlrtAddField(y, f_y, "fileTimeLo", i);
    c_u = r0->fileTimeHi;
    g_y = NULL;
    m3 = mxCreateNumericMatrix(1, 1, mxUINT32_CLASS, mxREAL);
    *(uint32_T *)mxGetData(m3) = c_u;
    emlrtAssign(&g_y, m3);
    emlrtAddField(y, g_y, "fileTimeHi", i);
    c_u = r0->mFileTimeLo;
    h_y = NULL;
    m3 = mxCreateNumericMatrix(1, 1, mxUINT32_CLASS, mxREAL);
    *(uint32_T *)mxGetData(m3) = c_u;
    emlrtAssign(&h_y, m3);
    emlrtAddField(y, h_y, "mFileTimeLo", i);
    c_u = r0->mFileTimeHi;
    i_y = NULL;
    m3 = mxCreateNumericMatrix(1, 1, mxUINT32_CLASS, mxREAL);
    *(uint32_T *)mxGetData(m3) = c_u;
    emlrtAssign(&i_y, m3);
    emlrtAddField(y, i_y, "mFileTimeHi", i);
  }

  emlrtAssign(&nameCaptureInfo, y);
  emlrtNameCapturePostProcessR2012a(emlrtAlias(nameCaptureInfo));
  return nameCaptureInfo;
}

void vibor_t_api(const mxArray * const prhs[20], const mxArray *plhs[3])
{
  emxArray_real_T *Pc;
  emxArray_real_T *Pg;
  emxArray_real_T *Pw;
  struct_T RC;
  emxArray_real_T *Sw1;
  emxArray_real_T *Sw2;
  emxArray_real_T *CL;
  emxArray_real_T *GL;
  emxArray_real_T *dVc;
  emxArray_real_T *dVg;
  emxArray_real_T *W1C;
  emxArray_real_T *W1G;
  emxArray_real_T *WoC;
  emxArray_real_T *WoG;
  real_T ndt;
  real_T fl;
  b_struct_T PR;
  real_T dt;
  real_T i;
  real_T j_ndt;
  boolean_T fl1;
  emlrtHeapReferenceStackEnterFcnR2012b(emlrtRootTLSGlobal);
  emxInit_real_T(&Pc, 1, &d_emlrtRTEI, TRUE);
  emxInit_real_T(&Pg, 1, &d_emlrtRTEI, TRUE);
  emxInit_real_T(&Pw, 1, &d_emlrtRTEI, TRUE);
  emxInitStruct_struct_T(&RC, &d_emlrtRTEI, TRUE);
  emxInit_real_T(&Sw1, 1, &d_emlrtRTEI, TRUE);
  emxInit_real_T(&Sw2, 1, &d_emlrtRTEI, TRUE);
  emxInit_real_T(&CL, 1, &d_emlrtRTEI, TRUE);
  b_emxInit_real_T(&GL, 2, &d_emlrtRTEI, TRUE);
  emxInit_real_T(&dVc, 1, &d_emlrtRTEI, TRUE);
  emxInit_real_T(&dVg, 1, &d_emlrtRTEI, TRUE);
  emxInit_real_T(&W1C, 1, &d_emlrtRTEI, TRUE);
  emxInit_real_T(&W1G, 1, &d_emlrtRTEI, TRUE);
  b_emxInit_real_T(&WoC, 2, &d_emlrtRTEI, TRUE);
  b_emxInit_real_T(&WoG, 2, &d_emlrtRTEI, TRUE);

  /* Marshall function inputs */
  ndt = emlrt_marshallIn(emlrtAliasP(prhs[0]), "ndt");
  fl = emlrt_marshallIn(emlrtAliasP(prhs[1]), "fl");
  c_emlrt_marshallIn(emlrtAlias(prhs[2]), "Pc", Pc);
  c_emlrt_marshallIn(emlrtAlias(prhs[3]), "Pg", Pg);
  c_emlrt_marshallIn(emlrtAlias(prhs[4]), "Pw", Pw);
  e_emlrt_marshallIn(emlrtAliasP(prhs[5]), "PR", &PR);
  j_emlrt_marshallIn(emlrtAliasP(prhs[6]), "RC", &RC);
  dt = emlrt_marshallIn(emlrtAliasP(prhs[7]), "dt");
  c_emlrt_marshallIn(emlrtAlias(prhs[8]), "Sw1", Sw1);
  c_emlrt_marshallIn(emlrtAlias(prhs[9]), "Sw2", Sw2);
  c_emlrt_marshallIn(emlrtAlias(prhs[10]), "CL", CL);
  o_emlrt_marshallIn(emlrtAlias(prhs[11]), "GL", GL);
  c_emlrt_marshallIn(emlrtAlias(prhs[12]), "dVc", dVc);
  c_emlrt_marshallIn(emlrtAlias(prhs[13]), "dVg", dVg);
  c_emlrt_marshallIn(emlrtAlias(prhs[14]), "W1C", W1C);
  c_emlrt_marshallIn(emlrtAlias(prhs[15]), "W1G", W1G);
  q_emlrt_marshallIn(emlrtAlias(prhs[16]), "WoC", WoC);
  q_emlrt_marshallIn(emlrtAlias(prhs[17]), "WoG", WoG);
  i = emlrt_marshallIn(emlrtAliasP(prhs[18]), "i");
  j_ndt = emlrt_marshallIn(emlrtAliasP(prhs[19]), "j_ndt");

  /* Invoke the target function */
  fl1 = vibor_t(&ndt, fl, Pc, Pg, Pw, &PR, &RC, dt, Sw1, Sw2, CL, GL, dVc, dVg,
                W1C, W1G, WoC, WoG, i, &j_ndt);

  /* Marshall function outputs */
  plhs[0] = emlrt_marshallOut(ndt);
  plhs[1] = emlrt_marshallOut(j_ndt);
  plhs[2] = b_emlrt_marshallOut(fl1);
  WoG->canFreeData = FALSE;
  emxFree_real_T(&WoG);
  WoC->canFreeData = FALSE;
  emxFree_real_T(&WoC);
  W1G->canFreeData = FALSE;
  emxFree_real_T(&W1G);
  W1C->canFreeData = FALSE;
  emxFree_real_T(&W1C);
  dVg->canFreeData = FALSE;
  emxFree_real_T(&dVg);
  dVc->canFreeData = FALSE;
  emxFree_real_T(&dVc);
  GL->canFreeData = FALSE;
  emxFree_real_T(&GL);
  CL->canFreeData = FALSE;
  emxFree_real_T(&CL);
  Sw2->canFreeData = FALSE;
  emxFree_real_T(&Sw2);
  Sw1->canFreeData = FALSE;
  emxFree_real_T(&Sw1);
  emxFreeStruct_struct_T(&RC);
  Pw->canFreeData = FALSE;
  emxFree_real_T(&Pw);
  Pg->canFreeData = FALSE;
  emxFree_real_T(&Pg);
  Pc->canFreeData = FALSE;
  emxFree_real_T(&Pc);
  emlrtHeapReferenceStackLeaveFcnR2012b(emlrtRootTLSGlobal);
}

/* End of code generation (vibor_t_api.c) */
