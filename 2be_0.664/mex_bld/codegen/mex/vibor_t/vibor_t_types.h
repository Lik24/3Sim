/*
 * vibor_t_types.h
 *
 * Code generation for function 'vibor_t'
 *
 * C source code generated on: Mon Jul 21 16:57:34 2014
 *
 */

#ifndef __VIBOR_T_TYPES_H__
#define __VIBOR_T_TYPES_H__

/* Include files */
#include "rtwtypes.h"

/* Type Definitions */
#ifndef typedef_ResolvedFunctionInfo
#define typedef_ResolvedFunctionInfo
typedef struct
{
    const char * context;
    const char * name;
    const char * dominantType;
    const char * resolved;
    uint32_T fileTimeLo;
    uint32_T fileTimeHi;
    uint32_T mFileTimeLo;
    uint32_T mFileTimeHi;
} ResolvedFunctionInfo;
#endif /*typedef_ResolvedFunctionInfo*/
#ifndef typedef_b_struct_T
#define typedef_b_struct_T
typedef struct
{
    real_T Nl;
    real_T as[3];
    real_T aw[3];
    real_T ts[3];
    real_T tw[3];
    real_T mu[4];
    real_T Ns;
    real_T Ta;
    real_T dt;
    real_T ndt;
    real_T Ro[5];
    real_T bet;
    real_T lam[5];
    real_T Cp[5];
    real_T dh;
    real_T Kc;
    real_T zc[5];
    real_T Bo[5];
    real_T Fc;
    real_T Fc2;
    real_T Sc2;
    real_T dl;
} b_struct_T;
#endif /*typedef_b_struct_T*/
#ifndef struct_emxArray__common
#define struct_emxArray__common
struct emxArray__common
{
    void *data;
    int32_T *size;
    int32_T allocatedSize;
    int32_T numDimensions;
    boolean_T canFreeData;
};
#endif /*struct_emxArray__common*/
#ifndef typedef_emxArray__common
#define typedef_emxArray__common
typedef struct emxArray__common emxArray__common;
#endif /*typedef_emxArray__common*/
#ifndef struct_emxArray_boolean_T
#define struct_emxArray_boolean_T
struct emxArray_boolean_T
{
    boolean_T *data;
    int32_T *size;
    int32_T allocatedSize;
    int32_T numDimensions;
    boolean_T canFreeData;
};
#endif /*struct_emxArray_boolean_T*/
#ifndef typedef_emxArray_boolean_T
#define typedef_emxArray_boolean_T
typedef struct emxArray_boolean_T emxArray_boolean_T;
#endif /*typedef_emxArray_boolean_T*/
#ifndef struct_emxArray_int32_T
#define struct_emxArray_int32_T
struct emxArray_int32_T
{
    int32_T *data;
    int32_T *size;
    int32_T allocatedSize;
    int32_T numDimensions;
    boolean_T canFreeData;
};
#endif /*struct_emxArray_int32_T*/
#ifndef typedef_emxArray_int32_T
#define typedef_emxArray_int32_T
typedef struct emxArray_int32_T emxArray_int32_T;
#endif /*typedef_emxArray_int32_T*/
#ifndef struct_emxArray_real_T
#define struct_emxArray_real_T
struct emxArray_real_T
{
    real_T *data;
    int32_T *size;
    int32_T allocatedSize;
    int32_T numDimensions;
    boolean_T canFreeData;
};
#endif /*struct_emxArray_real_T*/
#ifndef typedef_emxArray_real_T
#define typedef_emxArray_real_T
typedef struct emxArray_real_T emxArray_real_T;
#endif /*typedef_emxArray_real_T*/
#ifndef typedef_struct_T
#define typedef_struct_T
typedef struct
{
    emxArray_real_T *rc;
    emxArray_real_T *Arc;
    emxArray_real_T *Cr;
    emxArray_real_T *Cc;
    emxArray_real_T *Gr;
    emxArray_real_T *Gc;
    emxArray_real_T *ACr;
    emxArray_real_T *ACc;
    emxArray_real_T *AGr;
    emxArray_real_T *AGc;
    real_T na;
    real_T nc;
    real_T ng;
} struct_T;
#endif /*typedef_struct_T*/

#endif
/* End of code generation (vibor_t_types.h) */
