/*
 * tomo_model.i --
 *
 *-----------------------------------------------------------------------------
 */

if (is_func(plug_in)) plug_in, "tomo_model";
// setup_package, "tomo_model";

struct lf2d_parameters {
  double pixelSize;                        /* detector pixel size */
  double pixelStep;                        /* detector pixel step (>=
                                              pixelSize) */
  double voxelSize;                        /* voxel size */
  double sourceDistance;                   /* distance origin-source */
  double detectorDistance;                 /* distance origin-detector */
  double objectOffset1, objectOffset2;     /* offsets of object center
                                              w.r.t. to its geometrical
                                              center */
  long objectDimension1, objectDimension2; /* dimensions of object */
  long detectorDimension;                  /* size of detector */
}

extern getpid;

extern _lf2d_direct;
/* PROTOTYPE
   void lf2d_direct(pointer, double, double, double array, double array, int);
 */

func lf2d_direct(p, beta, offset, src, dst, clr)
/* DOCUMENT lf2d_direct, prm, beta, offset, src, dst, clr;

   SEE ALSO:
 */
{
  if (structof(p) != lf2d_parameters) error, "bad parameters";
  dims = dimsof(src);
  if (structof(src) != double || numberof(dims) != 3 ||
      dims(2) != p.objectDimension1 || dims(3) != p.objectDimension2) {
    error, "bad SRC";
  }
  dims = dimsof(dst);
  if (structof(dst) != double || numberof(dims) != 2 ||
      dims(2) != p.detectorDimension) {
    error, "bad DST";
  }
  _lf2d_direct, &p, beta, offset, src, dst, (clr?1n:0n);
  return dst;
}

extern _lf2d_line_integ_direct;
/* PROTOTYPE
   void lf2d_line_integ_direct(pointer, double, double, double array, double array, int);
 */

func lf2d_line_integ_direct(p, beta, offset, src, dst, clr)
/* DOCUMENT lf2d_direct, prm, beta, offset, src, dst, clr;

   SEE ALSO:
 */
{
  if (structof(p) != lf2d_parameters) error, "bad parameters";
  dims = dimsof(src);
  if (structof(src) != double || numberof(dims) != 3 ||
      dims(2) != p.objectDimension1 || dims(3) != p.objectDimension2) {
    error, "bad SRC";
  }
  dims = dimsof(dst);
  if (structof(dst) != double || numberof(dims) != 2 ||
      dims(2) != p.detectorDimension) {
    error, "bad DST";
  }
  _lf2d_line_integ_direct, &p, beta, offset, src, dst, (clr?1n:0n);
  return dst;
}

extern _lf2d_adjoint;
/* PROTOTYPE
   void lf2d_adjoint(pointer, double, double, double array, double array, int);
 */

func lf2d_adjoint(p, beta, offset, src, dst, clr)
/* DOCUMENT lf2d_adjoint, prm, beta, offset, src, dst, clr;

   SEE ALSO:
 */
{
  if (structof(p) != lf2d_parameters) error, "bad parameters";
  dims = dimsof(src);
  if (structof(src) != double || numberof(dims) != 2 ||
      dims(2) != p.detectorDimension) {
    error, "bad SRC";
  }
  dims = dimsof(dst);
  if (structof(dst) != double || numberof(dims) != 3 ||
      dims(2) != p.objectDimension1 || dims(3) != p.objectDimension2) {
    error, "bad DST";
  }
  _lf2d_adjoint, &p, beta, offset, src, dst, (clr?1n:0n);
  return dst;
}

extern _lf2d_line_integ_adjoint;
/* PROTOTYPE
   void lf2d_line_integ_adjoint(pointer, double, double, double array, double array, int);
 */

func lf2d_line_integ_adjoint(p, beta, offset, src, dst, clr)
/* DOCUMENT lf2d_line_integ_adjoint, prm, beta, offset, src, dst, clr;

   SEE ALSO:
 */
{
  if (structof(p) != lf2d_parameters) error, "bad parameters";
  dims = dimsof(src);
  if (structof(src) != double || numberof(dims) != 2 ||
      dims(2) != p.detectorDimension) {
    error, "bad SRC";
  }
  dims = dimsof(dst);
  if (structof(dst) != double || numberof(dims) != 3 ||
      dims(2) != p.objectDimension1 || dims(3) != p.objectDimension2) {
    error, "bad DST";
  }
  _lf2d_line_integ_adjoint, &p, beta, offset, src, dst, (clr?1n:0n);
  return dst;
}

extern _lf2d_sparser;
/* PROTOTYPE
   void lf2d_sparser(pointer, double, double, long, long, double array, int);
 */

func lf2d_sparser(p, beta, offset, i1, i2, dst, clr)
/* DOCUMENT lf2d_sparser, prm, beta, offset, i1, i2, dst, clr;

   SEE ALSO:
 */
{
  if (structof(p) != lf2d_parameters) error, "bad parameters";
  dims = dimsof(src);
  dims = dimsof(dst);
  if (structof(dst) != double || numberof(dims) != 2 ||
      dims(2) != p.detectorDimension) {
    error, "bad DST";
  }
  _lf2d_sparser, &p, beta, offset, i1, i2, dst, (clr?1n:0n);
  return dst;
}

extern _lf2d_line_integ_sparser;
/* PROTOTYPE
   void lf2d_line_integ_sparser(pointer, double, double, long, long, double array, int);
 */

func lf2d_line_integ_sparser(p, beta, offset, i1, i2, dst, clr)
/* DOCUMENT lf2d_line_integ_sparser, prm, beta, offset, i1, i2, dst, clr;

   SEE ALSO:
 */
{
  if (structof(p) != lf2d_parameters) error, "bad parameters";
  dims = dimsof(src);
  dims = dimsof(dst);
  if (structof(dst) != double || numberof(dims) != 2 ||
      dims(2) != p.detectorDimension) {
    error, "bad DST";
  }
  _lf2d_line_integ_sparser, &p, beta, offset, i1, i2, dst, (clr?1n:0n);
  return dst;
}

func lf2d_test1(beta, color=)
{
  if (is_void(beta)) beta = 0.0;
  objectDimension1 = 3;
  objectDimension2 = 4;
  detectorDimension = 400;
  p = lf2d_parameters(pixelSize=0.9e-3,
                      pixelStep=1.0e-3,
                      voxelSize=1.0e-2,
                      sourceDistance = 541e-3, /* distance origin-source */
                      detectorDistance = 408e-3, /* distance origin-detector */
                      objectOffset1 = 0.0,
                      objectOffset2 = 0.0,
                      objectDimension1 = objectDimension1,
                      objectDimension2 = objectDimension2,
                      detectorDimension = detectorDimension);

  voxel = array(double, objectDimension1, objectDimension2);
  voxel(objectDimension1/2, objectDimension2/2) = 1.1;
  pixel = array(double, detectorDimension);
  lf2d_direct, p, beta, 0.0, voxel, pixel;
  plg, pixel, indgen(1:detectorDimension), color=color;

  x = random_n(objectDimension1, objectDimension2);
  y = random_n(detectorDimension);
  Hx = array(double, detectorDimension);
  Hty = array(double, objectDimension1, objectDimension2);
  lf2d_direct, p, beta, 0.0, x, Hx;
  lf2d_adjoint, p, beta, 0.0, y, Hty;

  s1 = sum(y*Hx);
  s2 = sum(Hty*x);
  write, format = "<y,H.x>  = %g\n<H'.y,x> = %g\n", s1, s2;


}

func lf2d_new_operator(param, beta, offset)
{
  if (is_scalar(param) && structof(param) == lf2d_parameters) {
    param = param; // private copy
  } else {
    error, "bad parameters";
  }

  if (is_vector(beta) && identof(beta) <= Y_DOUBLE) {
    beta = double(beta(*)); // private copy
    numberOfAngles = numberof(beta);
  } else {
    error, "bad BETA";
  }

  if (is_void(offset)) {
    offset = array(0.0, numberOfAngles);
  } else if (is_vector(offset) && identof(offset) <= Y_DOUBLE
             && numberof(offset) == numberOfAngles) {
    offset = double(offset(*)); // private copy
  } else {
    error, "bad OFFSET";
  }

  ws = h_save(param, beta, offset, numberOfAngles,
              objectDim1 = param.objectDimension1,
              objectDim2 = param.objectDimension2,
              dataDim1 = param.detectorDimension,
              dataDim2 = numberOfAngles,
              nevals = 0);
  h_evaluator, ws, "_lf2d_apply";
  return ws;
}

func _lf2d_apply(ws, src, job)
{
  local beta; eq_nocopy, beta, ws.beta;
  local offset; eq_nocopy, offset, ws.offset;
  param = &ws.param;
  dataDim1 = ws.dataDim1;
  dataDim2 = ws.dataDim2; // also number of angles
  objectDim1 = ws.objectDim1;
  objectDim2 = ws.objectDim2;
  if (identof(src) > Y_DOUBLE) {
    error, "bad type for input array";
  }
  dims = dimsof(src);
  if (job) {
    /* Apply adjoint operator. */
    if (numberof(dims) != 3 ||
        dims(2) != dataDim1 || dims(3) != dataDim2) {
      error, "bad dimensions for input array";
    }
    dst = array(double, objectDim1, objectDim2);
    for (k = 1; k <= dataDim2; ++k) {
      _lf2d_adjoint, param, beta(k), offset(k), src(,k), dst, 0n;
    }
  } else {
    /* Apply direct operator. */
    if (numberof(dims) != 3 ||
        dims(2) != objectDim1 || dims(3) != objectDim2) {
      error, "bad dimensions for input array";
    }
    dst = array(double, dataDim1, dataDim2);
    tmp = array(double, dataDim1);
    for (k = 1; k <= dataDim2; ++k) {
      _lf2d_direct, param, beta(k), offset(k), src, tmp, 1n;
      dst(,k) = tmp;
    }
  }
  h_save, ws, nevals = (ws.nevals + 1);
  return dst;
}

func lf2d_new_cost_function(nil, param=, data=, weight=, beta=, offset=,
                            mu=)
{
  if (! is_void(nil)) error, "unexpected non-keyword argument";
  H = lf2d_new_operator(param, beta, offset);
  dims = dimsof(data);
  if (identof(data) > Y_DOUBLE) {
    error, "bad data type";
  }
  if (numberof(dims) != 3 || dims(2) != H.dataDim1 || dims(3) != H.dataDim2) {
    error, "bad data dimensions";
  }
  if (is_void(weight)) {
    weight = array(1.0, dims);
  } else {
    if (identof(weight) > Y_DOUBLE) {
      error, "bad weight type";
    }
    if (min(weight) < 0) {
      error, "bad weight value(s)";
    }
    if (max(weight) <= 0) {
      error, "no valid data";
    }
    dims = dimsof(weight);
    if (numberof(dims) != 3 || dims(2) != H.dataDim1 || dims(3) != H.dataDim2) {
      error, "bad weight dimensions";
    }
  }

  if (is_void(mu)) {
    mu = 0.0;
  } else if (is_scalar(mu) && identof(mu) <= Y_DOUBLE && mu >= 0.0) {
    mu = double (mu);
  } else {
    error, "bad regularization weight";
  }
  if (mu > 0.0 && ! is_func(rgl_totvar)) {
    include, "totvar.i", 1;
  }

  ws = h_new(H = H, data = double(data), weight = double(weight), mu = mu);
  h_evaluator, ws, "_lf2d_cost";
  return ws;
}

func _lf2d_cost(ws, x, &gx)
{
  r = ws.H(x) - ws.data;
  wr = ws.weight*r;
  fx = 0.5*sum(wr*r);
  gx = ws.H(wr, 1n);
  if (ws.mu > 0.0) {
    local g;
    mu = ws.mu;
    fx += mu*rgl_totvar(x, g, threshold=1e-6);
    gx = unref(gx) + mu*g;
  }
  return fx;
}

func lf2d_load_data(file)
/* DOCUMENT fdata = lf2d_load_data(file);
         or fdata = lf2d_load_data();

     Load 2D tomographic data (cone beam geometry) and return the
     corresponding cost function.

   SEE ALSO:
 */
{
  if (is_void(file)) file = "data/DataSet.yhd";
  db = yhd_restore(file);
  objectDims = dimsof(db.trueObject);
  dataDims = dimsof(db.noisyData);
  param = lf2d_parameters(pixelSize = db.pixelSize,
                          pixelStep = db.pixelStep,
                          voxelSize = db.voxelSize,
                          sourceDistance = db.sourceDistance,
                          detectorDistance = db.detectorDistance,
                          objectOffset1 = db.objectOffset(1),
                          objectOffset2 = db.objectOffset(2),
                          objectDimension1 = objectDims(2),
                          objectDimension2 = objectDims(3),
                          detectorDimension = dataDims(2));
  return lf2d_new_cost_function(param = param, data = db.noisyData,
                                beta = (is_void(db.beta)
                                        ? db.theta - pi/2
                                        : db.beta));
}

func lf2d_new_viewer(nil, win=, cmin=, cmax=, cmap=, cbar=, vert=,
                     format=)
{
  if (! is_void(nil)) error, "unexpected non-keyword argument";
  ws = h_new(win=win, cmin=cmin, cmax=cmax, cmap=cmap,
             cbar=cbar, vert=vert, format=format);
  h_evaluator, ws, "_lf2d_viewer";
  return ws;
}
func _lf2d_viewer(ws, x, ..)
{
  if (is_void(ws.win)) {
    win = current_window();
    if (win < 0) {
      for (win = 0; win < 64; ++win) {
        if (! window_exists(win)) {
          break;
        }
      }
      if (win >= 64) {
        win = 0;
      }
    }
    h_set, ws, win=win;
  }
  window, ws.win;
  if (! is_void(ws.cmap)) {
    cmap, ws.cmap;
  }
  fma;
  SCALE = 1e3;
  WATER_ABSORPTION_COEFFICIENT = 0.1928; // cm^{-1}
  tmp = (SCALE/WATER_ABSORPTION_COEFFICIENT)*x - SCALE;
  pl_img, tmp, cmin=ws.cmin, cmax=ws.cmax, cbar=ws.cbar,
    vert=ws.vert, format=ws.format;
  pause, 1;
}

func lf2d_reconst(costfunc, xguess, mu=, maxiter=, maxeval=, verb=,
                  xmin=, xmax=, ftol=, gtol=, mem=, viewer=)
{
  if (! is_func(rgl_totvar)) include, "totvar.i", 1;
  if (! is_func(op_mnb)) include, "OptimPack1.i", 1;
  if (is_void(xguess)) {
    x = array(double,
              costfunc.H.param.objectDimension1,
              costfunc.H.param.objectDimension2);
  } else {
    x = double(xguess);
  }
  h_set, costfunc, mu = (is_void(mu) ? 0.0 : double(mu));
  return op_mnb(costfunc, x, mem=mem, maxiter=maxiter, maxeval=maxeval,
                xmin=xmin, xmax=xmax, /* ftol=ftol, gtol=gtol,*/
                verb=verb, viewer=viewer);
}

func lf2d_test2
{
  cf = lf2d_load_data();
  x = lf2d_reconst(cf, mu=0.0001, maxiter=1000, maxeval=1500,
                   verb=1, xmin=, xmax=, ftol=, gtol=, mem=3,
                   viewer=lf2d_new_viewer(win=7,cmin=0,cmax=100,
                                          vert=1,cbar=1,cmap="Reds",
                                          format="%.0f"));
}

/*
 * Local Variables:
 * mode: Yorick
 * tab-width: 8
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * fill-column: 78
 * coding: utf-8
 * End:
 */
