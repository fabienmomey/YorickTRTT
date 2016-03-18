// include, "Ytrtt2D_dyn.i", 1;
mp_include, "Ytrtt2D_dyn.i";

// data_dir= "/home/momey/Recherche_Tomographie/Data/data_CLB_09-11-2011/img_1.3.46.423632.135428.1320854260.10/";
data_dir= "/home/momey/Data/data_CLB_09-11-2011/mvt2d/";
// data_dir= "/home/momey/Data/data_CLB_09-11-2011/mvt3d/";
DATA = yhd_restore(data_dir+"Data4TRTT");

ntheta = DATA.ndata;
Rsc = DATA.Rsc;
Rsd = DATA.Rsd;
theta = DATA.theta; //FIXME: the first five projections are static
ntheta = numberof(theta);
theta = theta(5:0);
ndata=numberof(theta);
local data; eq_nocopy, data, DATA.data(..,5:);


/* TOMOBJ */
nx = 380;
ny = 380;
nt = 22;
xoff = -45.0; //FIXME: in mm
yoff = -20.0; //FIXME: in mm
s_scl = 1.0; //FIXME: in mm
s_deg = 3;
size_footprint = 10000;

/* TOMDATA */
nv = DATA.nv;
voff = DATA.voff(5:0); //FIXME: in mm
uoff = DATA.uoff(5:0);
v_scl = DATA.v_scl; //FIXME: in mm
t_index = 0.0;

uoff = long(floor((uoff/v_scl) +0.5));

/* PROJECTOR */
mode = TRTT_FAN_BEAM;
type = TRTT_SPLINE_DRIVEN;

/* Noise features */                  
photon_flux = 1.e6; //FIXME: SNR² ~ photon_flux
read_noise = 5.;
var_noise_uniform = 1.;

/* Extract data */
data2d = array(double,nv,ndata);
midpanel = long(floor(nv/2));
for (k=1; k<=ndata; ++k) {
    data2d(,k) = data(midpanel-uoff(k),,k);
}

id1 = where(theta>pi);
id2 = where(theta<=pi);
theta(id2) = -1.0*theta(id2);
theta(id1) = 2*pi-1.0*theta(id1);

DATA=[];
data=[];

/* Dynamic parameters */
datesfile = open(data_dir+"dates");
dates = array(double,ntheta);
_n=read(datesfile, format="%f\n", dates);
close, datesfile;
dates /= 1000.;
dates = dates(5:0);

tacq = 116.; //FIXME: in s (see FRAME.DBF)
tcycl = 4.; //FIXME: in s
t_data_scl = 0.1825; //FIXME: (see FRAME.DBF) => mean acquisition time step
t_deg = 1;

Htomo = trtt_2D_create_whole_simu_system(data2d,
                                         nv,
                                         v_scl,
                                         nx,
                                         ny,
                                         nt,
                                         s_scl,
                                         tcycl,
                                         tacq,
                                         s_deg,
                                         t_deg,
                                         size_footprint,
                                         ndata,
                                         mode,
                                         type,
                                         Rsc,
                                         Rsd,
                                         dates=dates,
                                         xoff=xoff,
                                         yoff=yoff,
                                         voff=-voff, /*FIXME:FIXME: IL FAUT METTRE -VOFF !!!*/
                                         theta=theta,
                                         // var_noise_uniform=var_noise_uniform,
                                         // photon_flux=photon_flux,
                                         // read_noise=read_noise,
                                         // add_noise=add_noise,
                                         use_sparse_coefs=0n);

local Cops; eq_nocopy, Cops, Htomo.Cops;
Ck_list=Cops.Ck_list;

trtt_plot_sino, Htomo.Cops, 1, dyn=1n;
data = trtt_get_sino(Cops,dyn=1n);
 
/* GO RECONSTRUCTION ! */
dweights = array(double,nv,ndata);
dweights(1:511,)=1.0;
id_invalid = where(data == -1.0);
if (is_array(id_invalid))
    dweights(id_invalid) = 0.0;

eps1=1.0; eps2=1.e-6;
mu_s=1.e-6; mu_t=1.0;
// regulTV = h_new(weight=[mu_s, mu_s, mu_t], threshold = eps1, options = RGL_TOTVAR_ISOTROPIC);
regulTV = h_new(weight=[mu_s, mu_t], threshold = [eps1, eps2], flag_separable=1n);

XR = trtt_2D_optim_simu_launcher(Cops, trtt_cost_quadratic_mpy_opky, use_sparse_coefs=0n, x=array(double,nx,ny,nt), dweights=dweights, regulTV=regulTV, xmin=0.0, mem=5, viewer=1n, win_viewer=4, win_viewer2=60, maxiter=300, maxeval=300, verbose=1n);
trtt_plot_slice, XR.x, 16, 5;
for(k=1; k<=5*nt; ++k) {
    trtt_plot_slice, XR.x, (k%nt), 6; pause, 50;
}

// for (k=1; k<=nt; ++k) {mywindow, 3; fma; img_plot, h_get(Htomo,"XR_mus5p10-5_mut5p10-4_eps10-12").x(..,(k-1)%nt+1), cmin=0.0, cmax=4.2, cbar=1, vert=1, format="%.3f"; palette, "gray.gp"; pdf, swrite(format="dataCLB/22frames_spatio-temporal/slice%02d",k);}
