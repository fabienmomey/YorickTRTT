
#define SCOPE          extern
#define INTEGER        int
#define LUDEC_ROW      NAME(LU_dec_row)
#define LUDEC_COL      NAME(LU_dec_col)
#define LUSLV          NAME(LU_slv)
#define LUDET          NAME(LU_det)
#define LUINV          NAME(LU_inv)


#define STYLE       2 /* column-major */
#define REAL        double
#define NAME(base)  base##_c_d
#include "gauss.c"
#undef REAL
#undef STYLE
#undef NAME

#define STYLE       2 /* column-major */
#define REAL        float
#define NAME(base)  base##_c_f
#include "gauss.c"
#undef REAL
#undef STYLE
#undef NAME

#define STYLE       3 /* row-major */
#define REAL        double
#define NAME(base)  base##_r_d
#include "gauss.c"
#undef REAL
#undef STYLE
#undef NAME

#define STYLE       3 /* row-major */
#define REAL        float
#define NAME(base)  base##_r_f
#include "gauss.c"
#undef REAL
#undef STYLE
#undef NAME

