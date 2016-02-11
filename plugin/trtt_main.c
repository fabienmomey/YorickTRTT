#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <time.h>

#include "trtt_common.h"
#include "trtt_3D_projectors.h"

int main(int argc, char *argv[])
{
    double *xos = (double*)calloc(12, sizeof(double));
    double *xsd = (double*)calloc(12, sizeof(double));
    long nx = 256;
    long ny = 256;
    long nz = 256;
    double s_scl= 1.0;
    int s_deg= 3;
    long nu = 512;
    long nv = 512;
    double u_scl = 1.0;
    double v_scl = 1.0;
    int job = 0;
    double *in = (double*)calloc(16777216, sizeof(double));
    double *out = (double*)calloc(262144, sizeof(double));
    int res;
    
    xos[0] = 1.0; xos[5] = 1.0; xos[10] = 1.0;
    xsd[0] = 1.0; xsd[5] = 1.0; xsd[10] = 1.0;
    xsd[3] = 1.0;

    res = trtt_FB3D(out, in, nx, ny, nz, s_scl, s_deg, nu, nv, u_scl, v_scl, xos, xsd, job);

    free(in);
    free(out);
    free(xos);
    free(xsd);
}
