// include, "Ytrtt2D_dyn.i", 1;
mp_include, "Ytrtt2D_dyn.i";

/* TOMOBJ */
nx = 380;
ny = 380;
nt = 22;

Htomo=trtt_2D_load("dataCLB_new_14_01_2016/results_dataset10_380x380x22_pixel1mm");

local Cops; eq_nocopy, Cops, Htomo.Cops;
Ck_list=Cops.Ck_list;

nv=Cops.C0001.Ak.Yk.nv;
ndata=numberof(Ck_list);

// trtt_plot_sino, Htomo.Cops, 1, dyn=1n;
data = trtt_get_sino(Cops,dyn=1n);

/* GO RECONSTRUCTION ! */
dweights = array(double,nv,ndata);
dweights(1:511,)=1.0;
id_invalid = where(data == -1.0);
if (is_array(id_invalid))
    dweights(id_invalid) = 0.0;

// include, "png.i", 1;
// mask=img_read("mask.png");
// mask/=255.0;
// mask=1.0-mask;
// mask=mask(..,-:1:nt);
// mask*=0.05;

// eps1=1.; eps2=1.e-6;
// mu_s=1.e-6; mu_t=1.0;
// regulTV = h_new(weight=[mu_s, mu_s, mu_t], threshold = eps1, options = RGL_TOTVAR_ISOTROPIC);
// regulTV = h_new(weight=[mu_s, mu_t], threshold = [eps1, eps2], flag_separable=1n);
// XR = trtt_2D_optim_simu_launcher(Cops, trtt_cost_quadratic_mpy_opky, use_sparse_coefs=0n, x=array(double,nx,ny,nt), dweights=dweights, regulTV=regulTV, xmin=0.0, mem=3, viewer=1n, win_viewer=4, win_viewer2=60, maxiter=100, maxeval=100, verbose=1n);

eps1=1.e-6; eps2=1.e-6;
// eps=1.e-3;
eps_name="R_eps1_1e-6_eps2_1e-6";
if (is_void(h_get(Htomo,eps_name))) {
    h_set, Htomo, eps_name, h_new();
    trtt_2D_save, Htomo, "dataCLB_new_14_01_2016/results_dataset10_380x380x22_pixel1mm", overwrite=1n;
    Htomo=trtt_2D_load("dataCLB_new_14_01_2016/results_dataset10_380x380x22_pixel1mm");
    local Cops; eq_nocopy, Cops, Htomo.Cops;
    Ck_list=Cops.Ck_list;
}

name_mu_s=["1e2","5e1","2e1","1e1","5e0","2e0","1e0","5e-1","2e-1","1e-1","5e-2","1e-2","5e-3"];
name_mu_t=["5e2","1e2","5e1","2e1","1e1","5e0","2e0","1e0","5e-1","2e-1","1e-1"];
mu_s_s=[100.0,50.0,20.0,10.0,5.0,2.0,1.0,0.5,0.2,0.1,0.05,0.01,0.005];
mu_t_s=[500.0,100.0,50.0,20.0,10.0,5.0,2.0,1.0,0.5,0.2,0.1];
Nmu_s=numberof(mu_s_s);
Nmu_t=numberof(mu_t_s);

for (i=1;i<=Nmu_t;++i) {
    for (j=1;j<=Nmu_s;++j) {        
        mu_s=mu_s_s(j); mu_t=mu_t_s(i);
        if (is_void(h_get(h_get(Htomo,eps_name),swrite(format="XR_mus%s_mut%s",name_mu_s(j),name_mu_t(i))))) {
            write, format="%s\n", swrite(format="XR_mus%s_mut%s",name_mu_s(j),name_mu_t(i));
            regulTV = h_new(weight=[mu_s, mu_t], threshold = [eps1, eps2], flag_separable=1n);
            // regulTV = h_new(weight=[mu_s, mu_s, mu_t], threshold = eps, options = RGL_TOTVAR_ISOTROPIC);
            XR = trtt_2D_optim_simu_launcher(Cops, trtt_cost_quadratic_mpy_opky, use_sparse_coefs=0n, x=array(double,nx,ny,nt), dweights=dweights, regulTV=regulTV, xmin=0.0, mem=3, viewer=0n, win_viewer=4, win_viewer2=60, maxiter=2000, maxeval=2000, verbose=1n);
            h_set, h_get(Htomo,eps_name), swrite(format="XR_mus%s_mut%s",name_mu_s(j),name_mu_t(i)), XR;
            
            trtt_2D_save, Htomo, "dataCLB_new_14_01_2016/results_dataset10_380x380x22_pixel1mm", overwrite=1n;
            Htomo=trtt_2D_load("dataCLB_new_14_01_2016/results_dataset10_380x380x22_pixel1mm");
            local Cops; eq_nocopy, Cops, Htomo.Cops;
            Ck_list=Cops.Ck_list;
        }
    }
}

mp_exec, "if (!mp_rank) quit;";

// trtt_plot_slice, XR.x, 16, 5;
// for(k=1; k<=10*nt; ++k) {
//     trtt_plot_slice, XR.x, (k%nt), 6; pause, 100;
// }

// h_set, Htomo, "Reps1e-1", h_new();
// h_set, XR, n_iter=424;
// h_set, XR, n_eval=443;


// h_set, h_get(Htomo,"Reps1e-1"), "XR_mus15_mut15", XR;
// h_show, h_get(h_get(Htomo,"Reps1e-1"),"XR_mus10_mut10")
// h_set, h_get(h_get(Htomo,"Reps1e-1"),"XR_mus10_mut15"), n_iter=424;
// h_set, h_get(h_get(Htomo,"Reps1e-1"),"XR_mus10_mut15"), n_eval=443;

