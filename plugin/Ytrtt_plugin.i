TRTT_PLUGIN_PATH = TRTT_PATH+"plugin";

current_dirs=plug_dir(TRTT_PLUGIN_PATH);
if (is_func(plug_in)) plug_in, "Ytrtt";
    
/**********************************
 * COLUMN-MAJOR, DOUBLE PRECISION *
 **********************************/
extern _LU_dec_row_c_d;
/* PROTOTYPE
   int LU_dec_row_c_d(double array a, int n, int array pvt, double array ws);
*/
extern _LU_dec_col_c_d;
/* PROTOTYPE
   int LU_dec_col_c_d(double array a, int n, int array pvt, double array ws);
*/
extern _LU_slv_c_d;
/* PROTOTYPE
   int LU_slv_c_d(double array a, int n, int array pvt, double array b);
*/
extern _LU_inv_c_d;
/* PROTOTYPE
   int LU_inv_c_d(double array a, int n, int array pvt,
                  double array b, double array ws);
*/
extern _LU_det_c_d;
/* PROTOTYPE
   double LU_det_c_d(double array a, int n, int array pvt);
*/

/**********************************
 * COLUMN-MAJOR, SINGLE PRECISION *
 **********************************/
extern _LU_dec_row_c_f;
/* PROTOTYPE
   int LU_dec_row_c_f(double array a, int n, int array pvt, double array ws);
*/
extern _LU_dec_col_c_f;
/* PROTOTYPE
   int LU_dec_col_c_f(double array a, int n, int array pvt, double array ws);
*/
extern _LU_slv_c_f;
/* PROTOTYPE
   int LU_slv_c_f(double array a, int n, int array pvt, double array b);
*/
extern _LU_inv_c_f;
/* PROTOTYPE
   int LU_inv_c_f(double array a, int n, int array pvt,
                  double array b, double array ws);
*/
extern _LU_det_c_f;
/* PROTOTYPE
   double LU_det_c_f(double array a, int n, int array pvt);
*/

/*******************************
 * ROW-MAJOR, DOUBLE PRECISION *
 *******************************/
extern _LU_dec_row_r_d;
/* PROTOTYPE
   int LU_dec_row_r_d(double array a, int n, int array pvt, double array ws);
*/
extern _LU_dec_col_r_d;
/* PROTOTYPE
   int LU_dec_col_r_d(double array a, int n, int array pvt, double array ws);
*/
extern _LU_slv_r_d;
/* PROTOTYPE
   int LU_slv_r_d(double array a, int n, int array pvt, double array b);
*/
extern _LU_inv_r_d;
/* PROTOTYPE
   int LU_inv_r_d(double array a, int n, int array pvt,
                  double array b, double array ws);
*/
extern _LU_det_r_d;
/* PROTOTYPE
   double LU_det_r_d(double array a, int n, int array pvt);
*/

/*******************************
 * ROW-MAJOR, SINGLE PRECISION *
 *******************************/
extern _LU_dec_row_r_f;
/* PROTOTYPE
   int LU_dec_row_r_f(double array a, int n, int array pvt, double array ws);
*/
extern _LU_dec_col_r_f;
/* PROTOTYPE
   int LU_dec_col_r_f(double array a, int n, int array pvt, double array ws);
*/
extern _LU_slv_r_f;
/* PROTOTYPE
   int LU_slv_r_f(double array a, int n, int array pvt, double array b);
*/
extern _LU_inv_r_f;
/* PROTOTYPE
   int LU_inv_r_f(double array a, int n, int array pvt,
                  double array b, double array ws);
*/
extern _LU_det_r_f;
/* PROTOTYPE
   double LU_det_r_f(double array a, int n, int array pvt);
*/

/*---------------------------------------------------------------------------*/

ROW_MAJOR = 1;
COLUMN_MAJOR = 0;
func LU_dec_row(a,order)
{
  if ((ident = identof(a)) > Y_DOUBLE ||
      (dims = dimsof(a))(1) != 2 || (n = dims(2)) != dims(3)) {
    error, "expecting a real square matrix";
  }
  if (ident <= Y_FLOAT) {
    op = (order ? _LU_dec_row_r_f : _LU_dec_row_c_f);
    type = float;
  } else {
    op = (order ? _LU_dec_row_r_d : _LU_dec_row_c_d);
    type = double;
  }
  a = type(a);
  pvt = array(int, n + 1);
  status = op(a, n, pvt, array(type, n));
  if (status) error, "matrix is singular";
  return h_new(a=a, pvt=pvt, order=order, n=n, old=0, type=type);
}
func LU_dec_col(a,order)
{
  if ((ident = identof(a)) > Y_DOUBLE ||
      (dims = dimsof(a))(1) != 2 || (n = dims(2)) != dims(3)) {
    error, "expecting a real square matrix";
  }
  if (ident <= Y_FLOAT) {
    op = (order ? _LU_dec_col_r_f : _LU_dec_col_c_f);
    type = float;
  } else {
    op = (order ? _LU_dec_col_r_d : _LU_dec_col_c_d);
    type = double;
  }
  a = type(a);
  pvt = array(int, n + 1);
  status = op(a, n, pvt, array(type, n));
  if (status) error, "matrix is singular";
  return h_new(a=a, pvt=pvt, order=order, n=n, old=1, type=type);
}
func LU_slv(lu,b)
{
  if (numberof(b) != (n = lu.n)) error, "bad size for B";
  type = lu.type;
  order = lu.order;
  if (type == float) {
    op = (order ? _LU_slv_r_f : _LU_slv_c_f);
  } else {
    op = (order ? _LU_slv_r_d : _LU_slv_c_d);
  }
  b = type(b);
  status = op(lu.a, n, lu.pvt, b);
  if (status) error, "matrix is singular";
  return b;
}
func LU_det(lu)
{
  if (lu.type == float) {
    op = (lu.order ? _LU_det_r_f : _LU_det_c_f);
  } else {
    op = (lu.order ? _LU_det_r_d : _LU_det_c_d);
  }
  return op(lu.a, lu.n, lu.pvt);
}
func LU_inv(lu)
{
  if ((type = lu.type) == float) {
    op = (lu.order ? _LU_inv_r_f : _LU_inv_c_f);
  } else {
    op = (lu.order ? _LU_inv_r_d : _LU_inv_c_d);
  }
  n = lu.n;
  ainv = array(type, n, n);
  status = op(lu.a, n, lu.pvt, ainv, array(type, n));
  if (status) error, "matrix is singular";
  return ainv;
}

func LU_tests(n, niter=, order=)
{
  if (is_void(niter)) niter = 1;
  if (! order) order = 0;

  a = random(n,n) - 0.5;
  x = random(n) - 0.5;

  

  if (order == ROW_MAJOR) {
    b = a(+,)*x(+);
  } else {
    b = a(,+)*x(+);
  }
  
  t0 = t1 = array(double, 3);
  
  if (order == ROW_MAJOR) {
    at = transpose(a);
    timer, t0;
    for (k=1;k<=niter;++k){
      x0 = LUsolve(at,b);
    }
    timer, t1;
  } else {
    timer, t0;
    for (k=1;k<=niter;++k){
      x0 = LUsolve(a,b);
    }
    timer, t1;
  }
  write, format="LUsolve:             max(err) = %g  in %g ms\n",
    max(abs(x-x0)), (t1(3) - t0(3))*1E3/niter;
  
  //ainv = LUsolve(a);
  //s = SVdec(a);
  //adet = exp(sum(log(s))); // absolute value of the determinant
  timer, t0;
  for (k=1;k<=niter;++k){
    lu1 = LU_dec_row(a,order);
    x1 = LU_slv(lu1,b);
  }
  timer, t1;
  write, format="LU_dec_row + LU_slv: max(err) = %g  in %g ms\n",
    max(abs(x-x1)), (t1(3) - t0(3))*1E3/niter;
  
  timer, t0;
  for (k=1;k<=niter;++k){
    lu2 = LU_dec_col(a,order);
    x2 = LU_slv(lu2,b);
  }
  timer, t1;
  write, format="LU_dec_col + LU_slv: max(err) = %g  in %g ms\n",
    max(abs(x-x2)), (t1(3) - t0(3))*1E3/niter;

  if (order == ROW_MAJOR) {
    ainv0 = transpose(LUsolve(transpose(a)));
  } else {
    ainv0 = LUsolve(a);
  }
  ainv1 = LU_inv(lu1);
  ainv2 = LU_inv(lu2);
  max(abs(ainv1 - ainv0))/max(abs(ainv0));
  max(abs(ainv2 - ainv0))/max(abs(ainv0));
  if (n < 7) {
    pm, ainv0;
    pm, ainv1;
    pm, ainv2;
  }
}

/*---------------------------------------------------------------------------*/
/* TRTT_XFORM3D */

/* FIXME: the semantics should be: if DST is missing and called as a function,
   a new object is allocated, otherwise the operation is done "in-place".
*/

extern trtt_xform3d_new;
/* DOCUMENT obj = trtt_xform3d_new()

     This function creates a 3D coordinate transform.  The returned object
     can be indexed or used as function to apply the coordinate transform:

        obj()     // get the 12 coefficients of the transform
        obj(i)    // get the i-th coefficient of the transform (i is a scalar
                  // integer)

   FIXME: not yet implemented
        obj(v)    // apply the transform to vector v a 3-by-any array
                  // v(1,..) = X-coordinates, v(2,..) = Y-coordinates, etc.
        

   SEE ALSO: trtt_xform3d_move, trtt_xform3d_rotate, trtt_xform3d_scale,
             trtt_xform3d_invert, trtt_xform3d_apply, trtt_xform3d_combine.
 */

extern trtt_xform3d_move;
/* DOCUMENT trtt_xform3d_move, dst, tx, ty, tz, src;
         or trtt_xform3d_move, dst, tx, ty, tz;
         or trtt_xform3d_move, dst, t, src;
         or trtt_xform3d_move, dst, t;
         or dst = trtt_xform3d_move(tx, ty, tz, src);
         or dst = trtt_xform3d_move(t, src);

     This subroutine applies a translation to a 3D coordinate transform.  DST
     and SRC are the destination and source 3D coordinate transforms.  If SRC
     is omitted, the operation is done "in-place".  The translation can be
     specified by its 3 coordinates TX, TY, and TZ or by a 3-element vector T.
     When called as a function, DST is returned.

   SEE ALSO: trtt_xform3d_new
 */

extern trtt_xform3d_rotate;
/* DOCUMENT trtt_xform3d_rotate, dst, qw, qx, qy, qz, src;
         or trtt_xform3d_rotate, dst, qw, qx, qy, qz;
         or trtt_xform3d_rotate, dst, q, src;
         or trtt_xform3d_rotate, dst, q;
         or dst = trtt_xform3d_rotate(q, src);
         or dst = trtt_xform3d_rotate(qw, qx, qy, qz, src);
   
     This subroutine applies a rotation to a 3D coordinate transform.  DST and
     SRC are the destination and source 3D coordinate transforms.  If SRC is
     omitted, the operation is done "in-place".  QW, QX, QY, and QZ are the
     coordinates of the quaternion which specifies the rotation.  When called
     as a function, DST is returned.
 */

extern trtt_xform3d_invert;
/* DOCUMENT trtt_xform3d_invert, dst, src;
         or trtt_xform3d_invert, dst;
         or dst = trtt_xform3d_invert(src);
   
     This subroutine inverts a 3D coordinate transform.  DST and SRC are the
     destination and source 3D coordinate transforms.  If SRC is omitted, the
     operation is done "in-place".  When called as a function, DST is
     returned.
 */

extern trtt_xform3d_scale;
/* DOCUMENT trtt_xform3d_scale, dst, scl, src;
         or trtt_xform3d_scale, dst, scl;
         or dst = trtt_xform3d_scale(scl, src);
   
     This subroutine scales a 3D coordinate transform.  DST and SRC are the
     destination and source 3D coordinate transforms.  If SRC is omitted the
     operation is done "in-place".  SCL is the scaling factor.  When called as
     a function, DST is returned.
 */

extern trtt_xform3d_copy;
/* DOCUMENT trtt_xform3d_copy, dst, src;
         or trtt_xform3d_copy(src);
   
     This subroutine copies a 3D coordinate transform.  DST and SRC are the
     destination and source 3D coordinate transforms.  If SRC is omitted, a
     new coordinate transform object is created.  SCL is the scaling factor.
     When called as a function, DST is returned.
 */

extern trtt_xform3d_reset;
/* DOCUMENT trtt_xform3d_reset, obj;
   
     This subroutine resets the 3D coordinate transform OBJ to be the identity.
     When called as a function, OBJ is returned.
 */

extern trtt_xform3d_set_coefficients;
/* DOCUMENT trtt_xform3d_set_coefficients, obj, c;
   
     This subroutine set the coefficients of the 3D coordinate transform OBJ.
     C is an array of 12 values.  When called as a function, OBJ is returned.
 */

extern trtt_xform3d_combine;
/* DOCUMENT trtt_xform3d_combine, dst, left, right;
         or trtt_xform3d_combine(left, right);
   
     This subroutine combines the 3D coordinate transforms LEFT and RIGHT to
     produce a new 3D coordinate transform equivalent to applying the RIGHT
     transform.  The result is stored in DST (which can be LEFT and/or RIGHT
     to overwrite one of these transforms); if DST is not provided, a new 3D
     coordinate transform is created and returned.  When called as a function,
     the result transform is returned.
 */

func _trtt_xform3d_tests_message(msg, a1, a2, a3, a4)
{
  if (is_void(a1)) {
    return msg;
  } else if (is_void(a2)) {
    return swrite(format=msg, a1);
  } else if (is_void(a3)) {
    return swrite(format=msg, a1, a2);
  } else if (is_void(a4)) {
    return swrite(format=msg, a1, a2, a3);
  } else {
    return swrite(format=msg, a1, a2, a3, a4);
  }
}

func _trtt_xform3d_tests_failure(msg, a1, a2, a3, a4)
{
  extern nfailures;
  ++nfailures;
  write, format = "*** FAILURE: %s\n",
    _trtt_xform3d_tests_message(msg, a1, a2, a3, a4);
}

func _trtt_xform3d_tests_success(msg, a1, a2, a3, a4)
{
  write, format = "    success: %s\n",
    _trtt_xform3d_tests_message(msg, a1, a2, a3, a4);
}

func trtt_xform3d_tests
{
  TRUE = 1n;
  FALSE = 0n;
  nfailures = 0;
  fmt = "Test: %s\n";
  failure = _trtt_xform3d_tests_failure;
  success = _trtt_xform3d_tests_success;
  
  a = trtt_xform3d_new();
  b = trtt_xform3d_new();
  c = trtt_xform3d_new();
  identity = a();

  /* check set/get coefficients */
  write, format = fmt, "set/get coefficients";
  xa = 10.0*random_n(12);
  trtt_xform3d_set_coefficients, a, xa;
  for (i=-11; i<=12; ++i) {
    if (xa(i) != a(i)) {
      failure, "bad i-th coefficient value (i = %d)", i;
    }
  }
  
  /* check copy */
  write, format = fmt, "copy";
  b = trtt_xform3d_copy(a);
  if (max(abs(xa - a())) != 0.0) {
    failure, "contents of object A has changed";
  }
  if (max(abs(xa - b())) != 0.0) {
    failure, "objects B and A have different coefficients";
  }
  ref = trtt_xform3d_copy(c, a);
  if (max(abs(xa - c())) != 0.0) {
    failure, "objects C and A have different coefficients";
  }
  if (ref != c) {
    failure, "destination and returned objects are not the same";
  }

  /* check move */
  write, format = fmt, "move";
  t = 10.0*random_n(3);
  trtt_xform3d_move, b, t, a;
  if (max(abs(xa - a())) != 0.0) {
    failure, "contents of object A has changed";
  }
  trtt_xform3d_move, c, t(1), t(2), t(3), a;
  if (max(abs(xa - a())) != 0.0) {
    failure, "contents of object A has changed";
  }
  xb = b();  
  if (max(abs(xb - c())) != 0.0) {
    failure, "objects B and C have different coefficients";
  }
  // FIXME: check also in-place move

  /* check reset */
  write, format = fmt, "reset";
  ref = trtt_xform3d_reset(trtt_xform3d_copy(b, a));
  if (max(abs(identity - b())) != 0.0) {
    failure, "coefficients are not those of the identity";
  }
  if (max(abs(xa - a())) != 0.0) {
    failure, "contents of object A has changed";
  }
  if (ref != b) {
    failure, "destination and returned objects are not the same";
  }
  
  /* check that: invert(I) ~ I */
  write, format = fmt, "invert(I) ~ I";
  trtt_xform3d_reset, b;
  trtt_xform3d_invert, c, b;
  if (max(abs(identity - c())) != 0.0) {
    failure, "coefficients are not those of the identity";
  }
  
  
  /* check that: invert(move(A, t)) ~ move(OBJ, -t) */
  
  /* check that: invert(rotate(A, theta)) ~ rotate(A, -theta) */
  
  /* check that: invert(invert(A)) ~ A
   *             combine(A,invert(A)) ~ I
   *             combine(invert(A),A) ~ I
   */
  write, format = fmt, "invert(invert(A)) ~ A";
  trtt_xform3d_invert, b, a;
  trtt_xform3d_invert, c, b;
  err = max(abs(c() - a()));
  if (err > 1E-13) failure, "inversion error (errmax = %.3E)", err;
  write, format = fmt, "combine(A,invert(A)) ~ I";
  trtt_xform3d_combine, c, a, b;
  err = max(abs(identity - c()));
  if (err > 1E-14) {
    failure, "coefficients are not those of the identity (errmax = %.3E)", err;
  }
  write, format = fmt, "combine(invert(A),A) ~ I";
  trtt_xform3d_combine, c, b, a;
  err = max(abs(identity - c()));
  if (err > 1E-14) {
    failure, "coefficients are not those of the identity (errmax = %.3E)", err;
  }

  /* check that: scale(A,alpha) ~ alpha*A
   *             scale(scale(A,alpha),1/alpha) ~ A
   */
  write, format = fmt, "scale(A,alpha) ~ alpha*A";
  alpha = 3.111;
  trtt_xform3d_set_coefficients, a, xa;
  ref = trtt_xform3d_scale(b, alpha, a);
  if (max(abs(xa - a())) != 0.0) {
    failure, "contents of object A has changed";
  }
  if (ref != b) {
    failure, "destination and returned objects are not the same";
  }
  xb = b();
  err = max(abs(alpha*xa - xb));
  if (err > 0.0) {
    failure, "coefficients have not expected value (errmax = %.3E)", err;
  }
  write, format = fmt, "scale(scale(A,alpha),1/alpha) ~ A";
  ref = trtt_xform3d_scale(c, 1.0/alpha, b);
  if (max(abs(xb - b())) != 0.0) {
    failure, "contents of object B has changed";
  }
  if (ref != c) {
    failure, "destination and returned objects are not the same";
  }
  err = max(abs(xa - c()));
  if (err > 0.0) {
    failure, "coefficients have not expected value (errmax = %.3E)", err;
  }
  

    

  if (nfailures == 0) {
    summary = "OK: all tests successfully passed";
  } else {
    summary = swrite(format = "ERROR: %d failure(s)", nfailures);
  }
  write, format="*** %s\n", summary;
}

/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/* TRTT_PROJECTORS */

extern trtt_mvmult;
/* PROTOTYPE
int trtt_mvmult(double array, double array, long, double array, long array, long array)
*/

extern trtt_PB2D_cubic;
/* PROTOTYPE
int trtt_PB2D_cubic(double array, double array, long, long, double, int, long, double, double array, double array, int)
*/

extern trtt_FB2D_cubic;
/* PROTOTYPE
int trtt_FB2D_cubic(double array, double array, long, long, double, int, long, double, double array, double array, int)
*/

extern trtt_PB2D;
/* PROTOTYPE
int trtt_PB2D(double array, double array, long, long, double, int, double array, double array, long, double, double array, long, double, double array, double array, int)
*/

extern trtt_FB2D;
/* PROTOTYPE
int trtt_FB2D(double array, double array, long, long, double, int, double array, double array, long, double, double array, long, double, double array, double array, int)
*/

extern trtt_PB2D_sparser;
/* PROTOTYPE
int trtt_PB2D_sparser(double array, long, long, long, long, double, int, double array, double array, long, double, double array, long, double, double array, double array)
*/

extern trtt_FB2D_sparser;
/* PROTOTYPE
int trtt_FB2D_sparser(double array, long, long, long, long, double, int, double array, double array, long, double, double array, long, double, double array, double array)
*/

extern trtt_PB3D;
/* PROTOTYPE
int trtt_PB3D(double array, double array, long, long, long, double, int, double array, double array, long, double, double array, long, long, double, double, double array, double array, int)
*/

extern trtt_DD_PB3D;
/* PROTOTYPE
int trtt_DD_PB3D(double array, double array, long, long, long, double, long, long, double, double, double array, double array, int)
*/

extern trtt_CB3D;
/* PROTOTYPE
int trtt_CB3D(double array, double array, long, long, long, double, int, double array, double array, long, double, double array, long, long, double, double, double array, double array, int)
*/

extern trtt_PB3D_cubic;
/* PROTOTYPE
   int trtt_PB3D_cubic(double array, double array, long, long, long, double, int, long, long, double, double, double array, double array, int)
 */

extern trtt_CB3D_cubic;
/* PROTOTYPE
   int trtt_CB3D_cubic(double array, double array, long, long, long, double, int, long, long, double, double, double array, double array, int)
 */
