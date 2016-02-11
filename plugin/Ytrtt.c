#include <yapi.h>
#include "trtt_xform3d.h"

/* y_error is a no-return function.  This macro does explicit return
   to avoid compiler warnings.  Note that this macro is designed for
   built-in functions. */
#define ERROR(msg) do { y_error(msg); return; } while (0)

typedef struct _xform3d_instance xform3d_instance_t;
struct _xform3d_instance {
  trtt_xform3d_t *xform;
};

static void xform3d_free(void *);
static void xform3d_print(void *);
static void xform3d_eval(void *, int);
static void xform3d_extract(void *, char *);

static y_userobj_t xform3d_class = {
  "XForm3D",
  xform3d_free,
  xform3d_print,
  xform3d_eval,
  xform3d_extract,
  NULL
};

/* Get an instance of an XForm3D object at position iarg from the top of
   the stack. */
#define XFORM3D_GET(iarg) \
  ((xform3d_instance_t*)yget_obj(iarg,&xform3d_class))

/* Push a new XForm3D object on top of the stack. */
static xform3d_instance_t *xform3d_push(void)
{
  xform3d_instance_t *obj;
  obj = (xform3d_instance_t *)ypush_obj(&xform3d_class,
                                        sizeof(xform3d_instance_t));
  if ((obj->xform = trtt_xform3d_new()) == NULL) {
    yarg_drop(1);
    y_error("insufficent memory");
  }
  return obj;
}

static void xform3d_free(void *addr)
{
  xform3d_instance_t *obj = (xform3d_instance_t *)addr;
  if (obj->xform != NULL) {
    trtt_xform3d_delete(obj->xform);
  }
}

static void xform3d_print(void *addr)
{
  double c[12];
  xform3d_instance_t *obj = (xform3d_instance_t *)addr;
  int i;
  char buf[100];
  trtt_xform3d_get_coefficients(c, obj->xform);
  y_print("XForm3D:", 0);
  for (i = 0; i < 12; ++i) {
    if (i > 0 && (i%4) == 0) {
      y_print("        ", 0);
    }
    sprintf(buf, "  %12.3E", c[i]);
    y_print(buf, (i%4) == 3);
  }
}

static void xform3d_eval(void *addr, int argc)
{
  double c[12];
  xform3d_instance_t *obj = (xform3d_instance_t *)addr;
  long i;
  int typeid;
  if (argc == 1) {
    typeid = yarg_typeid(0);
    if (typeid == Y_VOID) {
      long dims[2];
      dims[0] = 1;
      dims[1] = 12;
      trtt_xform3d_get_coefficients(ypush_d(dims), obj->xform);
      return;
    } else if (yarg_rank(0) == 0 && typeid >= Y_CHAR && typeid <= Y_LONG) {
      /* Single argument is a scalar integer i, return i-th coefficients. */
      i = ygets_l(0);
      if (i <= 0) {
        i += 12;
      }
      if (i < 1 || i > 12) ERROR("out of range index");
      trtt_xform3d_get_coefficients(c, obj->xform);
      ypush_double(c[i - 1]);
      return;
    }
  }
  ERROR("expecting a single nil or integer argument");
}

static void xform3d_extract(void *addr, char *member)
{
  xform3d_instance_t *obj = (xform3d_instance_t *)addr;
  int c;
  long dims[2];

  if (member != NULL && (c = member[0]) != 0 && member[1] == 0) {
      switch (c) {
      case 'a':
	  dims[0] = 1;
	  dims[1] = 12;		     
	  trtt_xform3d_get_coefficients(ypush_d(dims), obj->xform);
	  return;
      }
  }
  ERROR("non-existing member");
}

/* this function implements the built-in routine "trtt_xform3d_new" */
void Y_trtt_xform3d_new(int argc)
{
  if (argc != 1 || ! yarg_nil(0)) ERROR("expecting a single nil argument");
  xform3d_push();
}

/* this function implements the built-in routine "trtt_xform3d_invert" */
void Y_trtt_xform3d_invert(int argc)
{
  xform3d_instance_t *src, *dst;
  if (argc == 1) {
    src = XFORM3D_GET(0);
    if (yarg_subroutine()) {
      /* operation will be performed in-place */
      dst = src;
    } else {
      /* create a new XForm3D object to store the result on top of the stack */
      dst = xform3d_push();
    }
  } else if (argc == 2) {
    src = XFORM3D_GET(0);
    dst = XFORM3D_GET(1);
  } else {
    ERROR("expecting one or two arguments");
  }
  trtt_xform3d_invert(dst->xform, src->xform);
  if (argc == 2) {
    /* left DST on top of the stack */
    yarg_drop(1);
  }
}

/* this function implements the built-in routine "xform3d_set_coefficients" */
void Y_trtt_xform3d_set_coefficients(int argc)
{
  long ntot;
  double *c;
  if (argc != 2) ERROR("expecting two arguments");
  c = ygeta_d(0, &ntot, NULL);
  if (ntot != 12) ERROR("expecting 12 coefficients");
  trtt_xform3d_set_coefficients(XFORM3D_GET(1)->xform, c);
  yarg_drop(1); /* left XForm3D object on top of the stack */
}

void Y_trtt_xform3d_move(int argc)
{
  double tx, ty, tz;
  const xform3d_instance_t *src;
  xform3d_instance_t *dst;
  long ntot;
  const double *t;
  int iarg = argc;
  int in_place = ((argc != 2 && argc != 4) || yarg_subroutine());
  if (argc == 2 || argc == 3) {
    if (in_place) { 
      dst = XFORM3D_GET(--iarg);
      t = ygeta_d(--iarg, &ntot, NULL);
      src = (iarg > 0 ? XFORM3D_GET(--iarg) : dst);
    } else {
      t = ygeta_d(--iarg, &ntot, NULL);
      src = XFORM3D_GET(--iarg);
      dst = xform3d_push();
    }
    if (ntot != 3) ERROR("expecting a vector of 3 elements");
    tx = t[0];
    ty = t[1];
    tz = t[2];
  } else if (argc == 4 || argc == 5) {
    if (in_place) { 
      dst = XFORM3D_GET(--iarg);
      tx = ygets_d(--iarg);
      ty = ygets_d(--iarg);
      tz = ygets_d(--iarg);
      src = (iarg > 0 ? XFORM3D_GET(--iarg) : dst);
    } else {
      tx = ygets_d(--iarg);
      ty = ygets_d(--iarg);
      tz = ygets_d(--iarg);
      src = XFORM3D_GET(--iarg);
      dst = xform3d_push();
    }
  } else {
    ERROR("bad number of arguments");
  }
  trtt_xform3d_move(dst->xform, tx, ty, tz, src->xform);
  if (in_place) yarg_drop(argc - 1); /* left DST on top of the stack */
}

void Y_trtt_xform3d_rotate(int argc)
{
  double qw, qx, qy, qz;
  const double *q;
  const xform3d_instance_t *src;
  xform3d_instance_t *dst;
  long ntot;
  int iarg = argc;
  int in_place = ((argc != 2 && argc != 5) || yarg_subroutine());
  if (argc == 2 || argc == 3) {
    if (in_place) {
      dst = XFORM3D_GET(--iarg);
      q = ygeta_d(--iarg, &ntot, NULL);
      src = (iarg > 0 ? XFORM3D_GET(--iarg) : dst);
    } else {
      q = ygeta_d(--iarg, &ntot, NULL);
      src = XFORM3D_GET(--iarg);
      dst = xform3d_push();
    }
    if (ntot != 4) ERROR("expecting a vector of 4 elements");
    qw = q[0];
    qx = q[1];
    qy = q[2];
    qz = q[3];
  } else if (argc == 5 || argc == 6) {
    if (in_place) {
      dst = XFORM3D_GET(--iarg);
      qw = ygets_d(--iarg);
      qx = ygets_d(--iarg);
      qy = ygets_d(--iarg);
      qz = ygets_d(--iarg);
      src = (iarg > 0 ? XFORM3D_GET(--iarg) : dst);
    } else {
      qw = ygets_d(--iarg);
      qx = ygets_d(--iarg);
      qy = ygets_d(--iarg);
      qz = ygets_d(--iarg);
      src = XFORM3D_GET(--iarg);
      dst = xform3d_push();
    }
  } else {
    ERROR("bad number of arguments");
  }
  trtt_xform3d_rotate(dst->xform, qw, qx, qy, qz, src->xform);
  if (in_place) yarg_drop(argc - 1); /* left DST on top of the stack */
}

void Y_trtt_xform3d_scale(int argc)
{
  const xform3d_instance_t *src;
  xform3d_instance_t *dst;
  double x_scl, y_scl, z_scl;
  int in_place, iarg = argc;
  if (argc != 4 && argc != 5) ERROR("bad number of arguments");
  in_place = (argc == 5 || yarg_subroutine());
  if (in_place) {
    dst = XFORM3D_GET(--iarg);
    x_scl = ygets_d(--iarg);
    y_scl = ygets_d(--iarg);
    z_scl = ygets_d(--iarg);
    src = (iarg > 0 ? XFORM3D_GET(--iarg) : dst);
  } else {
    src = XFORM3D_GET(--iarg);
    x_scl = ygets_d(--iarg);
    y_scl = ygets_d(--iarg);
    z_scl = ygets_d(--iarg);
    dst = xform3d_push();
  }
  trtt_xform3d_scale(dst->xform, x_scl, y_scl, z_scl, src->xform);
  if (in_place) yarg_drop(argc - 1); /* left DST on top of the stack */
}

void Y_trtt_xform3d_reset(int argc)
{
  if (argc != 1) ERROR("expecting a single argument");
  trtt_xform3d_reset(XFORM3D_GET(0)->xform);
}

void Y_trtt_xform3d_copy(int argc)
{
  const xform3d_instance_t *src;
  xform3d_instance_t *dst;
  if (argc == 1) {
    src = XFORM3D_GET(0);
    dst = xform3d_push();
  } else if (argc == 2) {
    src = XFORM3D_GET(0);
    dst = XFORM3D_GET(1);
  } else {
    ERROR("bad number of arguments");
  }
  trtt_xform3d_copy(dst->xform, src->xform);
  if (argc > 1) yarg_drop(argc - 1); /* left DST on top of the stack */
}

void Y_trtt_xform3d_combine(int argc)
{
  xform3d_instance_t *dst;
  const xform3d_instance_t *left, *right;
  if (argc != 2 && argc != 3) ERROR("bad number of arguments");
  right = XFORM3D_GET(0);
  left = XFORM3D_GET(1);
  dst = (argc == 3 ? XFORM3D_GET(2) : xform3d_push());
  trtt_xform3d_combine(dst->xform, left->xform, right->xform);
  if (argc == 3) yarg_drop(2); /* left DST on top of the stack */
}

