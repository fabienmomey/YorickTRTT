// include, "Ytrtt2D_dyn.i", 1;
mp_include, "Ytrtt2D_dyn.i";

data_dir= "/home/momey/Recherche_Tomographie/Data/data_CLB_09-11-2011/img_1.3.46.423632.135428.1320854260.10/";
// data_dir= "/home/momey/Data/data_CLB_09-11-2011/mvt2d/";
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
data = trtt_get_sino(Cops,dyn=1n);
 
/* GO RECONSTRUCTION ! */
dweights = array(double,nv,ndata);
dweights(1:511,)=1.0;
id_invalid = where(data == -1.0);
if (is_array(id_invalid))
    dweights(id_invalid) = 0.0;

cost_mat=trtt_2D_optim_define_cost_mat(Cops, dweights);

count=10;
/***** BENCHMARK DATA-FIDELITY WITHOUT MPY *****/
write, "/***********************************************/";
write, "/***** BENCHMARK DATA-FIDELITY WITHOUT MPY *****/";
x=random(nx,ny,nt);
dims = dimsof(x);
nt = dims(0);
Cops = cost_mat.Cops;
WCyT = cost_mat.WCyT;
Ck_list = h_get(Cops,"Ck_list");
ndata = numberof(Ck_list);
timer_start;
for (i=1;i<=count;i++) {
    fx = 0.0;
    gx = 0.0;
    for (k=1; k<=ndata; ++k) {
        Ck = h_get(Cops, Ck_list(k));
        Ak = Ck.Ak;
        Bk = Ck.Bk;
        Yk = Ak.Yk;
        yk = Yk.pixels;
        Wk = Yk.Wk;    
        Cxk = Ck(x);
        Cxmyk = Cxk-yk;
        WCxmyk = Wk(Cxmyk);
        fx += 0.5*sum(Cxmyk*WCxmyk);   
        if (!exist_sparse_coefs) {
            gx += Ck(Wk(Cxk),1)(*);
        } else {
            MCk = Ak.MCk;
            gx += Bk(MCk(Wk(Cxk),1),1)(*);
        }
    }
    gx -= WCyT(*);
    gx = reform(gx,dims);
}
timer_elapsed(count);
/***********************************************/
/***** BENCHMARK DATA-FIDELITY WITH MPY *****/
write, "/********************************************/";
write, "/***** BENCHMARK DATA-FIDELITY WITH MPY *****/";
x=random(nx,ny,nt);
timer_start;
for (i=1;i<=count;i++) {
    gx = 0.0;
    fx=trtt_2D_mpy_dyn_eval(x, gx, Cops);
}
timer_elapsed(count);
/********************************************/
/******* BENCHMARK REGULARIZATION SEP *******/
write, "/********************************************/";
write, "/******* BENCHMARK REGULARIZATION SEP *******/";
eps1=1.0;
eps2=1.e-6;
mu_s =1.0;
mu_t =1.0;
regulTV = h_new(weight=[mu_s, mu_t], threshold = [eps1,eps2], flag_separable=1n);
x=random(nx,ny,nz,nt);
gx=array(double,nx,ny,nz,nt);
timer_start;
for (i=1;i<=count;i++) {
    fx = rgl_mixed_ndpt(mu_s, eps1, mu_t, eps2, x, gx, 1n);
}
timer_elapsed(count);
/********************************************/
/****** BENCHMARK REGULARIZATION GLOB *******/
write, "/********************************************/";
write, "/****** BENCHMARK REGULARIZATION GLOB *******/";
eps2=1.e-6;
mu_s =1.0;
mu_t =1.0;
regulTV = h_new(weight=[mu_s, mu_s, mu_t], threshold = eps2);
x=random(nx,ny,nt);
gx=array(double,nx,ny,nt);
timer_start;
for (i=1;i<=count;i++) {
    fx = rgl_totvar(x, gx, weight=[mu_s,mu_s,mu_t], threshold=eps2, options=RGL_TOTVAR_ISOTROPIC);
}
timer_elapsed(count);
/********************************************/

// /* MA MACHINE */
// /***********************************************/
// /***** BENCHMARK DATA-FIDELITY WITHOUT MPY *****/
// [28.4482,0.0116027,28.4659]
// /********************************************/
// /***** BENCHMARK DATA-FIDELITY WITH MPY *****/
// [11.1094,0.0084014,11.1182]
// /********************************************/
// /******* BENCHMARK REGULARIZATION SEP *******/
// [0.0585501,0,0.0585374]
// /********************************************/
// /****** BENCHMARK REGULARIZATION GLOB *******/
// [0.0543726,0,0.054364]
