// include, "Ytrtt4D.i", 1;
mp_include, "Ytrtt4D.i";

// data_dir= "/home/momey/Recherche_Tomographie/Data/data_CLB_patient_02-07-2012/img_1.3.46.423632.141000.1169042526.68/";
// data_dir= "/home/momey/Data/data_CLB_09-11-2011/mvt2d/";
// data_dir= "/home/momey/Data/data_CLB_09-11-2011/mvt3d/";
data_dir= "/home/momey/Data/data_CLB_patient_02-07-2012/";
DATA = yhd_restore(data_dir+"Data4TRTT");

ntheta = DATA.ndata;
Rsc = DATA.Rsc;
Rsd = DATA.Rsd;
theta = DATA.theta; //FIXME: the first five projections are static
ntheta = numberof(theta);
theta = theta(5:-1);
ndata=numberof(theta);
local data; eq_nocopy, data, DATA.data(..,5:-1);

/* TOMOBJ */
nx = 120;
ny = 120;
nz = 72;
nt = 13;
x_off = 72.0; //FIXME: in mm
y_off = 0.0; //FIXME: in mm
z_off = 0.0; //FIXME: in mm
s_scl = 4.0; //FIXME: in mm
s_deg = 2;
size_footprint = 10000;

/* TOMDATA */
nv = DATA.nv;
nu = DATA.nu;
u_off = DATA.uoff(5:-1);
v_off = DATA.voff(5:-1); //FIXME: in mm
u_scl = DATA.u_scl; //FIXME: in mm
v_scl = DATA.v_scl; //FIXME: in mm
t_index = 0.0;

/* PROJECTOR */
mode = TRTT_CONE_BEAM;
type = TRTT_SPLINE_DRIVEN;

/* Noise features */                  
photon_flux = 1.e6; //FIXME: SNR² ~ photon_flux
read_noise = 5.;
var_noise_uniform = 1.;

/* Extract data */
// data3d = array(double,nu,nv,ndata);
// for (k=1; k<=ndata; ++k) {
//     data3d(..,k) = transpose(data(..,k))(::-1,);
// }

id1 = where(theta>pi);
id2 = where(theta<=pi);
if (is_array(id2)) theta(id2) = -1.0*theta(id2);
if (is_array(id1)) theta(id1) = 2*pi-1.0*theta(id1);

DATA=[];

/***********************************************************************/
/**** JEU DE DONNEES FANTÔME MÉCANIQUE *********************************/
/***********************************************************************/
// data=data(146:367,..); //FIXME: Troncature détecteur en u (filtre anti-diffusion)
// nu = dimsof(data)(2);

// datesfile = open(data_dir+"dates");
// dates = array(double,ntheta);
// _n=read(datesfile, format="%f\n", dates);
// close, datesfile;
// dates /= 1000.;
// dates = dates(5:0);

// tacq = 116.; //FIXME: in s (see FRAME.DBF)
// tcycl = 4.; //FIXME: in s
// t_data_scl = 0.1825; //FIXME: (see FRAME.DBF) => mean acquisition time step
// t_deg = 1;

/***********************************************************************/
/**** JEU DE DONNEES PATIENT : MANIPULATION DU SIGNAL RESPIRATOIRE *****/
/***********************************************************************/

/* Dynamic parameters */
datesfile = open(data_dir+"dates");
dates = array(double,ntheta);
_n=read(datesfile, format="%f\n", dates);
close, datesfile;
dates /= 1000.;
dates = dates(5:-1);

//FIXME: période moyenne ~ 2.4 secondes
sigfile = open(data_dir+"phase.sig");
sig = array(double,ntheta);
_n=read(sigfile, format="%f\n", sig);
sig=sig(5:-1);

/* sur-échantillonnage du signal temporel */
ndata2=1000*ndata;
dates2=span(dates(1),dates(0),ndata2);
wd = spl_interp_coefs_new_spline((indgen(0:ndata2-1)/(ndata2-1.))*(ndata-1.), ndata, 3);
sig2=wd(spl_spline_samples_to_coefficients(sig,3));
avgsig=avg(sig2);
sig2=sig2-avgsig;
sig=sig-avgsig;

/* recherche des minimums et maximums locaux */
mins=[];
maxs=[];
idmins=[];
idmaxs=[];

for(k=2;k<=ndata2-1;++k) {
    if (sig2(k-1)>sig2(k) & sig2(k+1)>=sig2(k)) {
        grow, mins, sig2(k);
        grow, idmins, k;
    }
    else if (sig2(k-1)<sig2(k) & sig2(k+1)<=sig2(k)) {
        grow, maxs, sig2(k);
        grow, idmaxs, k;
    }
}

// window, 1; fma; plg, sig, dates; // plg, sig2, dates2, color="red";
// window, 1; plg, zeros, dates2(indices), type="none", marks=1, marker='\4', color="blue";

/* Période moyenne */
tcycl = 0.5*(avg(abs((dates2(idmins))(dif)))+avg(abs((dates2(idmaxs))(dif))));
half_tcycl = 0.5*tcycl;
dates3=dates2(idmins);
dates4=dates2(idmaxs);
grow, dates4, dates3;
dates4=dates4(sort(dates4));

/* tronquer les périodes incomplètes */
idtronc = where(dates>=dates3(1) & dates<dates3(0));
if(is_array(idtronc)) {
    dates = dates(idtronc);
    theta = theta(idtronc);
    data = data(..,idtronc);
    u_off = u_off(idtronc);
    v_off = v_off(idtronc);
    ndata = numberof(theta);
    sig3=sig(idtronc);
}

/* Normalisation des dates par les demi-périodes moyennes */
nextr = numberof(dates4);
newdates = array(double,ndata);
t_old=dates4(1);
for (u=1; u<=nextr-1; ++u) {
    Ttemp = dates4(u+1)-dates4(u);
    idTtemp = where(dates>=dates4(u) & dates<dates4(u+1));
    if(is_array(idTtemp)) {
        newdates(idTtemp)=(dates(idTtemp)-t_old)*(half_tcycl/Ttemp)+(u-1)*half_tcycl;
    }
    t_old = dates4(u+1);
}

// window, 1; plg, sig3, dates3(1)+newdates, color="blue";

/* temporal spline degree */
t_deg = 1;

/***********************************************************************/
/***********************************************************************/
/***********************************************************************/

Htomo = trtt_4D_create_whole_simu_system(data,
                                         nu,
                                         nv,
                                         u_scl,
                                         v_scl,
                                         nx,
                                         ny,
                                         nz,
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
                                         dates=newdates,
                                         x_off=x_off,
                                         y_off=y_off,
                                         z_off=z_off,
                                         u_off=-u_off,  /*FIXME:FIXME: IL FAUT METTRE -U_OFF !!!*/
                                         v_off=-v_off,  /*FIXME:FIXME: IL FAUT METTRE -V_OFF !!!*/
                                         theta=theta);
                                         // var_noise_uniform=var_noise_uniform,
                                         // photon_flux=photon_flux,
                                         // read_noise=read_noise,
                                         // add_noise=add_noise);

local Cops; eq_nocopy, Cops, Htomo.Cops;
data = trtt3D_get_sino(Cops,dyn=1n);

dweights = array(1.0,nu,nv,ndata);
id_invalid = where(data < 0.0);
if (is_array(id_invalid))
    dweights(id_invalid) = 0.0;
dweights(,1:5,)=0.0;
dweights(,nv-4:nv,)=0.0;
dweights(1:5,..)=0.0;
dweights(nu-4:nu,..)=0.0;

cost_mat=trtt_4D_optim_define_cost_mat(Cops, dweights);

count=10;
/***** BENCHMARK DATA-FIDELITY WITHOUT MPY *****/
write, "/***********************************************/";
write, "/***** BENCHMARK DATA-FIDELITY WITHOUT MPY *****/";
x=random(nx,ny,nz,nt);
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
x=random(nx,ny,nz,nt);
timer_start;
for (i=1;i<=count;i++) {
    gx = 0.0;
    fx=trtt_4D_mpy_eval(x, gx, Cops);
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
eps1=1.0;
eps2=1.e-6;
mu_s =1.0;
mu_t =1.0;
regulTV = h_new(weight=[mu_s, mu_s, mu_t], threshold =eps2, flag_separable=1n);
x=random(nx,ny,nz,nt);
gx=array(double,nx,ny,nz,nt);
timer_start;
for (i=1;i<=count;i++) {
    fx = rgl_totvar(x, gx, weight=[mu_s,mu_s,mu_s,mu_t], threshold=eps2, options=RGL_TOTVAR_ISOTROPIC);
}
timer_elapsed(count);
/********************************************/

mp_exec, "if (!mp_rank) quit;";

// /* KUBILAI 12 CPUs */
// /***********************************************/
// /***** BENCHMARK DATA-FIDELITY WITHOUT MPY *****/
// [1198.77,33.8312,1231.66]
// /********************************************/
// /***** BENCHMARK DATA-FIDELITY WITH MPY *****/
// [151.68,3.8404,155.414]
// /********************************************/
// /******* BENCHMARK REGULARIZATION SEP *******/
// [0.5668,0,0.566145]
// /********************************************/
// /****** BENCHMARK REGULARIZATION GLOB *******/
// [0.9284,0.024,0.951959]

// /* MA MACHINE ZBOOK CORE i7 4 CPUs */
// /***********************************************/
// /***** BENCHMARK DATA-FIDELITY WITHOUT MPY *****/
// [981.324,10.9862,992.175]
// /********************************************/
// /***** BENCHMARK DATA-FIDELITY WITH MPY *****/
// [339.31,2.75615,342.02]
// /********************************************/
// /******* BENCHMARK REGULARIZATION SEP *******/
// [0.556449,0.0011999,0.557735]
// /********************************************/
// /****** BENCHMARK REGULARIZATION GLOB *******/
// [0.813479,0.0063985,0.82013]

