/* codger-generated yorick package wrapper file */
#include "play.h"
#include "ydata.h"

/*----------------begin tomo_model.i */
extern BuiltIn Y_getpid;
extern BuiltIn Y__lf2d_direct;

extern void lf2d_direct(void *, double , double , double *, double *, 
  int );
void
Y__lf2d_direct(int n)
{
  if (n!=6) YError("_lf2d_direct takes exactly 6 arguments");
  lf2d_direct(yarg_sp(5), yarg_sd(4), yarg_sd(3), yarg_d(2,0), 
    yarg_d(1,0), yarg_si(0));
}

extern BuiltIn Y__lf2d_line_integ_direct;

extern void lf2d_line_integ_direct(void *, double , double , 
  double *, double *, int );
void
Y__lf2d_line_integ_direct(int n)
{
  if (n!=6) YError("_lf2d_line_integ_direct takes exactly 6 arguments");
  lf2d_line_integ_direct(yarg_sp(5), yarg_sd(4), yarg_sd(3), 
    yarg_d(2,0), yarg_d(1,0), yarg_si(0));
}

extern BuiltIn Y__lf2d_adjoint;

extern void lf2d_adjoint(void *, double , double , double *, 
  double *, int );
void
Y__lf2d_adjoint(int n)
{
  if (n!=6) YError("_lf2d_adjoint takes exactly 6 arguments");
  lf2d_adjoint(yarg_sp(5), yarg_sd(4), yarg_sd(3), yarg_d(2,0), 
    yarg_d(1,0), yarg_si(0));
}

extern BuiltIn Y__lf2d_line_integ_adjoint;

extern void lf2d_line_integ_adjoint(void *, double , double , 
  double *, double *, int );
void
Y__lf2d_line_integ_adjoint(int n)
{
  if (n!=6) YError("_lf2d_line_integ_adjoint takes exactly 6 arguments");
  lf2d_line_integ_adjoint(yarg_sp(5), yarg_sd(4), yarg_sd(3), 
    yarg_d(2,0), yarg_d(1,0), yarg_si(0));
}

extern BuiltIn Y__lf2d_sparser;

extern void lf2d_sparser(void *, double , double , long , long , 
  double *, int );
void
Y__lf2d_sparser(int n)
{
  if (n!=7) YError("_lf2d_sparser takes exactly 7 arguments");
  lf2d_sparser(yarg_sp(6), yarg_sd(5), yarg_sd(4), yarg_sl(3), 
    yarg_sl(2), yarg_d(1,0), yarg_si(0));
}

extern BuiltIn Y__lf2d_line_integ_sparser;

extern void lf2d_line_integ_sparser(void *, double , double , 
  long , long , double *, int );
void
Y__lf2d_line_integ_sparser(int n)
{
  if (n!=7) YError("_lf2d_line_integ_sparser takes exactly 7 arguments");
  lf2d_line_integ_sparser(yarg_sp(6), yarg_sd(5), yarg_sd(4), 
    yarg_sl(3), yarg_sl(2), yarg_d(1,0), yarg_si(0));
}


/*----------------list include files */

static char *y0_includes[] = {
  "tomo_model.i",
  0,
  0
};

/*----------------collect pointers and names */

static BuiltIn *y0_routines[] = {
  &Y_getpid,
  &Y__lf2d_direct,
  &Y__lf2d_line_integ_direct,
  &Y__lf2d_adjoint,
  &Y__lf2d_line_integ_adjoint,
  &Y__lf2d_sparser,
  &Y__lf2d_line_integ_sparser,
  0
};

static void *y0_values[] = {
  0
};

static char *y0_names[] = {
  "getpid",
  "_lf2d_direct",
  "_lf2d_line_integ_direct",
  "_lf2d_adjoint",
  "_lf2d_line_integ_adjoint",
  "_lf2d_sparser",
  "_lf2d_line_integ_sparser",
  0
};

/*----------------define package initialization function */

PLUG_EXPORT char *yk_tomo_model(char ***,
                         BuiltIn ***, void ***, char ***);
static char *y0_pkgname = "tomo_model";

char *
yk_tomo_model(char ***ifiles,
       BuiltIn ***code, void ***data, char ***varname)
{
  *ifiles = y0_includes;
  *code = y0_routines;
  *data = y0_values;
  *varname = y0_names;
  return y0_pkgname;
}
