/* codger-generated yorick package wrapper file */
#include "play.h"
#include "ydata.h"

/*----------------begin Ytrtt_plugin.i */
extern BuiltIn Y__LU_dec_row_c_d;

extern int LU_dec_row_c_d(double *, int , int *, double *);
void
Y__LU_dec_row_c_d(int n)
{
  if (n!=4) YError("_LU_dec_row_c_d takes exactly 4 arguments");
  PushIntValue(LU_dec_row_c_d(yarg_d(3,0), yarg_si(2), yarg_i(1,0), 
    yarg_d(0,0)));
}

extern BuiltIn Y__LU_dec_col_c_d;

extern int LU_dec_col_c_d(double *, int , int *, double *);
void
Y__LU_dec_col_c_d(int n)
{
  if (n!=4) YError("_LU_dec_col_c_d takes exactly 4 arguments");
  PushIntValue(LU_dec_col_c_d(yarg_d(3,0), yarg_si(2), yarg_i(1,0), 
    yarg_d(0,0)));
}

extern BuiltIn Y__LU_slv_c_d;

extern int LU_slv_c_d(double *, int , int *, double *);
void
Y__LU_slv_c_d(int n)
{
  if (n!=4) YError("_LU_slv_c_d takes exactly 4 arguments");
  PushIntValue(LU_slv_c_d(yarg_d(3,0), yarg_si(2), yarg_i(1,0), 
    yarg_d(0,0)));
}

extern BuiltIn Y__LU_inv_c_d;

extern int LU_inv_c_d(double *, int , int *, double *, double *);
void
Y__LU_inv_c_d(int n)
{
  if (n!=5) YError("_LU_inv_c_d takes exactly 5 arguments");
  PushIntValue(LU_inv_c_d(yarg_d(4,0), yarg_si(3), yarg_i(2,0), 
    yarg_d(1,0), yarg_d(0,0)));
}

extern BuiltIn Y__LU_det_c_d;

extern double LU_det_c_d(double *, int , int *);
void
Y__LU_det_c_d(int n)
{
  if (n!=3) YError("_LU_det_c_d takes exactly 3 arguments");
  PushDoubleValue(LU_det_c_d(yarg_d(2,0), yarg_si(1), yarg_i(0,0)));
}

extern BuiltIn Y__LU_dec_row_c_f;

extern int LU_dec_row_c_f(double *, int , int *, double *);
void
Y__LU_dec_row_c_f(int n)
{
  if (n!=4) YError("_LU_dec_row_c_f takes exactly 4 arguments");
  PushIntValue(LU_dec_row_c_f(yarg_d(3,0), yarg_si(2), yarg_i(1,0), 
    yarg_d(0,0)));
}

extern BuiltIn Y__LU_dec_col_c_f;

extern int LU_dec_col_c_f(double *, int , int *, double *);
void
Y__LU_dec_col_c_f(int n)
{
  if (n!=4) YError("_LU_dec_col_c_f takes exactly 4 arguments");
  PushIntValue(LU_dec_col_c_f(yarg_d(3,0), yarg_si(2), yarg_i(1,0), 
    yarg_d(0,0)));
}

extern BuiltIn Y__LU_slv_c_f;

extern int LU_slv_c_f(double *, int , int *, double *);
void
Y__LU_slv_c_f(int n)
{
  if (n!=4) YError("_LU_slv_c_f takes exactly 4 arguments");
  PushIntValue(LU_slv_c_f(yarg_d(3,0), yarg_si(2), yarg_i(1,0), 
    yarg_d(0,0)));
}

extern BuiltIn Y__LU_inv_c_f;

extern int LU_inv_c_f(double *, int , int *, double *, double *);
void
Y__LU_inv_c_f(int n)
{
  if (n!=5) YError("_LU_inv_c_f takes exactly 5 arguments");
  PushIntValue(LU_inv_c_f(yarg_d(4,0), yarg_si(3), yarg_i(2,0), 
    yarg_d(1,0), yarg_d(0,0)));
}

extern BuiltIn Y__LU_det_c_f;

extern double LU_det_c_f(double *, int , int *);
void
Y__LU_det_c_f(int n)
{
  if (n!=3) YError("_LU_det_c_f takes exactly 3 arguments");
  PushDoubleValue(LU_det_c_f(yarg_d(2,0), yarg_si(1), yarg_i(0,0)));
}

extern BuiltIn Y__LU_dec_row_r_d;

extern int LU_dec_row_r_d(double *, int , int *, double *);
void
Y__LU_dec_row_r_d(int n)
{
  if (n!=4) YError("_LU_dec_row_r_d takes exactly 4 arguments");
  PushIntValue(LU_dec_row_r_d(yarg_d(3,0), yarg_si(2), yarg_i(1,0), 
    yarg_d(0,0)));
}

extern BuiltIn Y__LU_dec_col_r_d;

extern int LU_dec_col_r_d(double *, int , int *, double *);
void
Y__LU_dec_col_r_d(int n)
{
  if (n!=4) YError("_LU_dec_col_r_d takes exactly 4 arguments");
  PushIntValue(LU_dec_col_r_d(yarg_d(3,0), yarg_si(2), yarg_i(1,0), 
    yarg_d(0,0)));
}

extern BuiltIn Y__LU_slv_r_d;

extern int LU_slv_r_d(double *, int , int *, double *);
void
Y__LU_slv_r_d(int n)
{
  if (n!=4) YError("_LU_slv_r_d takes exactly 4 arguments");
  PushIntValue(LU_slv_r_d(yarg_d(3,0), yarg_si(2), yarg_i(1,0), 
    yarg_d(0,0)));
}

extern BuiltIn Y__LU_inv_r_d;

extern int LU_inv_r_d(double *, int , int *, double *, double *);
void
Y__LU_inv_r_d(int n)
{
  if (n!=5) YError("_LU_inv_r_d takes exactly 5 arguments");
  PushIntValue(LU_inv_r_d(yarg_d(4,0), yarg_si(3), yarg_i(2,0), 
    yarg_d(1,0), yarg_d(0,0)));
}

extern BuiltIn Y__LU_det_r_d;

extern double LU_det_r_d(double *, int , int *);
void
Y__LU_det_r_d(int n)
{
  if (n!=3) YError("_LU_det_r_d takes exactly 3 arguments");
  PushDoubleValue(LU_det_r_d(yarg_d(2,0), yarg_si(1), yarg_i(0,0)));
}

extern BuiltIn Y__LU_dec_row_r_f;

extern int LU_dec_row_r_f(double *, int , int *, double *);
void
Y__LU_dec_row_r_f(int n)
{
  if (n!=4) YError("_LU_dec_row_r_f takes exactly 4 arguments");
  PushIntValue(LU_dec_row_r_f(yarg_d(3,0), yarg_si(2), yarg_i(1,0), 
    yarg_d(0,0)));
}

extern BuiltIn Y__LU_dec_col_r_f;

extern int LU_dec_col_r_f(double *, int , int *, double *);
void
Y__LU_dec_col_r_f(int n)
{
  if (n!=4) YError("_LU_dec_col_r_f takes exactly 4 arguments");
  PushIntValue(LU_dec_col_r_f(yarg_d(3,0), yarg_si(2), yarg_i(1,0), 
    yarg_d(0,0)));
}

extern BuiltIn Y__LU_slv_r_f;

extern int LU_slv_r_f(double *, int , int *, double *);
void
Y__LU_slv_r_f(int n)
{
  if (n!=4) YError("_LU_slv_r_f takes exactly 4 arguments");
  PushIntValue(LU_slv_r_f(yarg_d(3,0), yarg_si(2), yarg_i(1,0), 
    yarg_d(0,0)));
}

extern BuiltIn Y__LU_inv_r_f;

extern int LU_inv_r_f(double *, int , int *, double *, double *);
void
Y__LU_inv_r_f(int n)
{
  if (n!=5) YError("_LU_inv_r_f takes exactly 5 arguments");
  PushIntValue(LU_inv_r_f(yarg_d(4,0), yarg_si(3), yarg_i(2,0), 
    yarg_d(1,0), yarg_d(0,0)));
}

extern BuiltIn Y__LU_det_r_f;

extern double LU_det_r_f(double *, int , int *);
void
Y__LU_det_r_f(int n)
{
  if (n!=3) YError("_LU_det_r_f takes exactly 3 arguments");
  PushDoubleValue(LU_det_r_f(yarg_d(2,0), yarg_si(1), yarg_i(0,0)));
}

extern BuiltIn Y_trtt_xform3d_new;
extern BuiltIn Y_trtt_xform3d_move;
extern BuiltIn Y_trtt_xform3d_rotate;
extern BuiltIn Y_trtt_xform3d_invert;
extern BuiltIn Y_trtt_xform3d_scale;
extern BuiltIn Y_trtt_xform3d_copy;
extern BuiltIn Y_trtt_xform3d_reset;
extern BuiltIn Y_trtt_xform3d_set_coefficients;
extern BuiltIn Y_trtt_xform3d_combine;
extern BuiltIn Y_trtt_mvmult;

extern int trtt_mvmult(double *, double *, long , double *, long *, 
  long *);
void
Y_trtt_mvmult(int n)
{
  if (n!=6) YError("trtt_mvmult takes exactly 6 arguments");
  PushIntValue(trtt_mvmult(yarg_d(5,0), yarg_d(4,0), yarg_sl(3), 
    yarg_d(2,0), yarg_l(1,0), yarg_l(0,0)));
}

extern BuiltIn Y_trtt_PB2D_cubic;

extern int trtt_PB2D_cubic(double *, double *, long , long , 
  double , int , long , double , double *, double *, int );
void
Y_trtt_PB2D_cubic(int n)
{
  if (n!=11) YError("trtt_PB2D_cubic takes exactly 11 arguments");
  PushIntValue(trtt_PB2D_cubic(yarg_d(10,0), yarg_d(9,0), yarg_sl(8), 
    yarg_sl(7), yarg_sd(6), yarg_si(5), yarg_sl(4), yarg_sd(3), 
    yarg_d(2,0), yarg_d(1,0), yarg_si(0)));
}

extern BuiltIn Y_trtt_FB2D_cubic;

extern int trtt_FB2D_cubic(double *, double *, long , long , 
  double , int , long , double , double *, double *, int );
void
Y_trtt_FB2D_cubic(int n)
{
  if (n!=11) YError("trtt_FB2D_cubic takes exactly 11 arguments");
  PushIntValue(trtt_FB2D_cubic(yarg_d(10,0), yarg_d(9,0), yarg_sl(8), 
    yarg_sl(7), yarg_sd(6), yarg_si(5), yarg_sl(4), yarg_sd(3), 
    yarg_d(2,0), yarg_d(1,0), yarg_si(0)));
}

extern BuiltIn Y_trtt_PB2D;

extern int trtt_PB2D(double *, double *, long , long , double , 
  int , double *, double *, long , double , double *, long , 
  double , double *, double *, int );
void
Y_trtt_PB2D(int n)
{
  if (n!=16) YError("trtt_PB2D takes exactly 16 arguments");
  PushIntValue(trtt_PB2D(yarg_d(15,0), yarg_d(14,0), yarg_sl(13), yarg_sl(12), 
    yarg_sd(11), yarg_si(10), yarg_d(9,0), yarg_d(8,0), yarg_sl(7), 
    yarg_sd(6), yarg_d(5,0), yarg_sl(4), yarg_sd(3), 
    yarg_d(2,0), yarg_d(1,0), yarg_si(0)));
}

extern BuiltIn Y_trtt_FB2D;

extern int trtt_FB2D(double *, double *, long , long , double , 
  int , double *, double *, long , double , double *, long , 
  double , double *, double *, int );
void
Y_trtt_FB2D(int n)
{
  if (n!=16) YError("trtt_FB2D takes exactly 16 arguments");
  PushIntValue(trtt_FB2D(yarg_d(15,0), yarg_d(14,0), yarg_sl(13), yarg_sl(12), 
    yarg_sd(11), yarg_si(10), yarg_d(9,0), yarg_d(8,0), yarg_sl(7), 
    yarg_sd(6), yarg_d(5,0), yarg_sl(4), yarg_sd(3), 
    yarg_d(2,0), yarg_d(1,0), yarg_si(0)));
}

extern BuiltIn Y_trtt_PB2D_sparser;

extern int trtt_PB2D_sparser(double *, long , long , long , long , 
  double , int , double *, double *, long , double , double *, 
  long , double , double *, double *);
void
Y_trtt_PB2D_sparser(int n)
{
  if (n!=16) YError("trtt_PB2D_sparser takes exactly 16 arguments");
  PushIntValue(trtt_PB2D_sparser(yarg_d(15,0), yarg_sl(14), yarg_sl(13), 
    yarg_sl(12), yarg_sl(11), yarg_sd(10), yarg_si(9), yarg_d(8,0), 
    yarg_d(7,0), yarg_sl(6), yarg_sd(5), yarg_d(4,0), yarg_sl(3), 
    yarg_sd(2), yarg_d(1,0), yarg_d(0,0)));
}

extern BuiltIn Y_trtt_FB2D_sparser;

extern int trtt_FB2D_sparser(double *, long , long , long , long , 
  double , int , double *, double *, long , double , double *, 
  long , double , double *, double *);
void
Y_trtt_FB2D_sparser(int n)
{
  if (n!=16) YError("trtt_FB2D_sparser takes exactly 16 arguments");
  PushIntValue(trtt_FB2D_sparser(yarg_d(15,0), yarg_sl(14), yarg_sl(13), 
    yarg_sl(12), yarg_sl(11), yarg_sd(10), yarg_si(9), yarg_d(8,0), 
    yarg_d(7,0), yarg_sl(6), yarg_sd(5), yarg_d(4,0), yarg_sl(3), 
    yarg_sd(2), yarg_d(1,0), yarg_d(0,0)));
}

extern BuiltIn Y_trtt_PB3D;

extern int trtt_PB3D(double *, double *, long , long , long , 
  double , int , double *, double *, long , double , double *, 
  long , long , double , double , double *, double *, int );
void
Y_trtt_PB3D(int n)
{
  if (n!=19) YError("trtt_PB3D takes exactly 19 arguments");
  PushIntValue(trtt_PB3D(yarg_d(18,0), yarg_d(17,0), yarg_sl(16), yarg_sl(15), 
    yarg_sl(14), yarg_sd(13), yarg_si(12), yarg_d(11,0), yarg_d(10,0), 
    yarg_sl(9), yarg_sd(8), yarg_d(7,0), yarg_sl(6), 
    yarg_sl(5), yarg_sd(4), yarg_sd(3), yarg_d(2,0), 
    yarg_d(1,0), yarg_si(0)));
}

extern BuiltIn Y_trtt_DD_PB3D;

extern int trtt_DD_PB3D(double *, double *, long , long , long , 
  double , long , long , double , double , double *, double *, 
  int );
void
Y_trtt_DD_PB3D(int n)
{
  if (n!=13) YError("trtt_DD_PB3D takes exactly 13 arguments");
  PushIntValue(trtt_DD_PB3D(yarg_d(12,0), yarg_d(11,0), yarg_sl(10), 
    yarg_sl(9), yarg_sl(8), yarg_sd(7), yarg_sl(6), yarg_sl(5), 
    yarg_sd(4), yarg_sd(3), yarg_d(2,0), yarg_d(1,0), yarg_si(0)));
}

extern BuiltIn Y_trtt_CB3D;

extern int trtt_CB3D(double *, double *, long , long , long , 
  double , int , double *, double *, long , double , double *, 
  long , long , double , double , double *, double *, int );
void
Y_trtt_CB3D(int n)
{
  if (n!=19) YError("trtt_CB3D takes exactly 19 arguments");
  PushIntValue(trtt_CB3D(yarg_d(18,0), yarg_d(17,0), yarg_sl(16), yarg_sl(15), 
    yarg_sl(14), yarg_sd(13), yarg_si(12), yarg_d(11,0), yarg_d(10,0), 
    yarg_sl(9), yarg_sd(8), yarg_d(7,0), yarg_sl(6), 
    yarg_sl(5), yarg_sd(4), yarg_sd(3), yarg_d(2,0), 
    yarg_d(1,0), yarg_si(0)));
}

extern BuiltIn Y_trtt_PB3D_cubic;

extern int trtt_PB3D_cubic(double *, double *, long , long , 
  long , double , int , long , long , double , double , double *, 
  double *, int );
void
Y_trtt_PB3D_cubic(int n)
{
  if (n!=14) YError("trtt_PB3D_cubic takes exactly 14 arguments");
  PushIntValue(trtt_PB3D_cubic(yarg_d(13,0), yarg_d(12,0), yarg_sl(11), 
    yarg_sl(10), yarg_sl(9), yarg_sd(8), yarg_si(7), yarg_sl(6), 
    yarg_sl(5), yarg_sd(4), yarg_sd(3), yarg_d(2,0), yarg_d(1,0), 
    yarg_si(0)));
}

extern BuiltIn Y_trtt_CB3D_cubic;

extern int trtt_CB3D_cubic(double *, double *, long , long , 
  long , double , int , long , long , double , double , double *, 
  double *, int );
void
Y_trtt_CB3D_cubic(int n)
{
  if (n!=14) YError("trtt_CB3D_cubic takes exactly 14 arguments");
  PushIntValue(trtt_CB3D_cubic(yarg_d(13,0), yarg_d(12,0), yarg_sl(11), 
    yarg_sl(10), yarg_sl(9), yarg_sd(8), yarg_si(7), yarg_sl(6), 
    yarg_sl(5), yarg_sd(4), yarg_sd(3), yarg_d(2,0), yarg_d(1,0), 
    yarg_si(0)));
}


/*----------------list include files */

static char *y0_includes[] = {
  "Ytrtt_plugin.i",
  0,
  0
};

/*----------------collect pointers and names */

static BuiltIn *y0_routines[] = {
  &Y__LU_dec_row_c_d,
  &Y__LU_dec_col_c_d,
  &Y__LU_slv_c_d,
  &Y__LU_inv_c_d,
  &Y__LU_det_c_d,
  &Y__LU_dec_row_c_f,
  &Y__LU_dec_col_c_f,
  &Y__LU_slv_c_f,
  &Y__LU_inv_c_f,
  &Y__LU_det_c_f,
  &Y__LU_dec_row_r_d,
  &Y__LU_dec_col_r_d,
  &Y__LU_slv_r_d,
  &Y__LU_inv_r_d,
  &Y__LU_det_r_d,
  &Y__LU_dec_row_r_f,
  &Y__LU_dec_col_r_f,
  &Y__LU_slv_r_f,
  &Y__LU_inv_r_f,
  &Y__LU_det_r_f,
  &Y_trtt_xform3d_new,
  &Y_trtt_xform3d_move,
  &Y_trtt_xform3d_rotate,
  &Y_trtt_xform3d_invert,
  &Y_trtt_xform3d_scale,
  &Y_trtt_xform3d_copy,
  &Y_trtt_xform3d_reset,
  &Y_trtt_xform3d_set_coefficients,
  &Y_trtt_xform3d_combine,
  &Y_trtt_mvmult,
  &Y_trtt_PB2D_cubic,
  &Y_trtt_FB2D_cubic,
  &Y_trtt_PB2D,
  &Y_trtt_FB2D,
  &Y_trtt_PB2D_sparser,
  &Y_trtt_FB2D_sparser,
  &Y_trtt_PB3D,
  &Y_trtt_DD_PB3D,
  &Y_trtt_CB3D,
  &Y_trtt_PB3D_cubic,
  &Y_trtt_CB3D_cubic,
  0
};

static void *y0_values[] = {
  0
};

static char *y0_names[] = {
  "_LU_dec_row_c_d",
  "_LU_dec_col_c_d",
  "_LU_slv_c_d",
  "_LU_inv_c_d",
  "_LU_det_c_d",
  "_LU_dec_row_c_f",
  "_LU_dec_col_c_f",
  "_LU_slv_c_f",
  "_LU_inv_c_f",
  "_LU_det_c_f",
  "_LU_dec_row_r_d",
  "_LU_dec_col_r_d",
  "_LU_slv_r_d",
  "_LU_inv_r_d",
  "_LU_det_r_d",
  "_LU_dec_row_r_f",
  "_LU_dec_col_r_f",
  "_LU_slv_r_f",
  "_LU_inv_r_f",
  "_LU_det_r_f",
  "trtt_xform3d_new",
  "trtt_xform3d_move",
  "trtt_xform3d_rotate",
  "trtt_xform3d_invert",
  "trtt_xform3d_scale",
  "trtt_xform3d_copy",
  "trtt_xform3d_reset",
  "trtt_xform3d_set_coefficients",
  "trtt_xform3d_combine",
  "trtt_mvmult",
  "trtt_PB2D_cubic",
  "trtt_FB2D_cubic",
  "trtt_PB2D",
  "trtt_FB2D",
  "trtt_PB2D_sparser",
  "trtt_FB2D_sparser",
  "trtt_PB3D",
  "trtt_DD_PB3D",
  "trtt_CB3D",
  "trtt_PB3D_cubic",
  "trtt_CB3D_cubic",
  0
};

/*----------------define package initialization function */

PLUG_EXPORT char *yk_Ytrtt(char ***,
                         BuiltIn ***, void ***, char ***);
static char *y0_pkgname = "Ytrtt";

char *
yk_Ytrtt(char ***ifiles,
       BuiltIn ***code, void ***data, char ***varname)
{
  *ifiles = y0_includes;
  *code = y0_routines;
  *data = y0_values;
  *varname = y0_names;
  return y0_pkgname;
}
