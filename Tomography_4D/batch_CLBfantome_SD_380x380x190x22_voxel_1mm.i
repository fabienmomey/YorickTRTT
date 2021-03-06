mp_include, "Ytrtt4D.i";
// include, "Ytrtt4D.i", 1;

Htomo_name =  "./results/CLBfantome_SD_380x380x190x22_voxel_1mm";
Htomo_reconst_name = "./results/CLBfantome_SD_380x380x190x22_voxel_1mm.rec";
mu_s = 1.0;
mu_t = 1.0;
eps = 1.e-2; // eps1=1.0; eps2=1.e-6;
XRname = "XRglob_mus1e0_mut1e0_eps1e-2"; // XRname = "XRsepa_mus1e0_mut1e0_eps1_1e0_eps2_1e-6";

if (!is_void(open(Htomo_name, "rb", 1))) {
    write, format="---------- %s reloaded ----------\n", Htomo_name;
    Htomo = trtt_4D_load(Htomo_name);
    Htomo_reconst = yhd_restore(Htomo_reconst_name);
    local dweights; eq_nocopy, dweights, Htomo.Cops.dweights;
} else {
    write, format="---------- %s created ----------\n", Htomo_name;
    /* PATIENT */
    // data_dir= "/home/momey/Data/data_CLB_patient_02-07-2012/";
    /* FANT�ME M�CANIQUE */
    data_dir= "/home/momey/Data/data_CLB_09-11-2011/mvt3d/";
    
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
    nx = 380;
    ny = 380;
    nz = 190;
    nt = 22;
    x_off = -45.0; //FIXME: in mm
    y_off = -20.0; //FIXME: in mm
    z_off = 0.0; //FIXME: in mm
    s_scl = 1.0; //FIXME: in mm
    s_deg = 3;
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
    photon_flux = 1.e6; //FIXME: SNR� ~ photon_flux
    read_noise = 5.;
    var_noise_uniform = 1.;
    
    id1 = where(theta>pi);
    id2 = where(theta<=pi);
    if (is_array(id2)) theta(id2) = -1.0*theta(id2);
    if (is_array(id1)) theta(id1) = 2*pi-1.0*theta(id1);
    
    DATA=[];
    
    /***********************************************************************/
    /**** JEU DE DONNEES FANT�ME M�CANIQUE *********************************/
    /***********************************************************************/
    data=data(146:367,..); //FIXME: Troncature d�tecteur en u (filtre anti-diffusion)
    nu = dimsof(data)(2);
    
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
    
    /***********************************************************************/
    /**** JEU DE DONNEES PATIENT : MANIPULATION DU SIGNAL RESPIRATOIRE *****/
    /***********************************************************************/
    
    // /* Dynamic parameters */
    // datesfile = open(data_dir+"dates");
    // dates = array(double,ntheta);
    // _n=read(datesfile, format="%f\n", dates);
    // close, datesfile;
    // dates /= 1000.;
    // dates = dates(5:-1);
    
    // //FIXME: p�riode moyenne ~ 2.4 secondes
    // sigfile = open(data_dir+"phase.sig");
    // sig = array(double,ntheta);
    // _n=read(sigfile, format="%f\n", sig);
    // sig=sig(5:-1);
    
    // /* sur-�chantillonnage du signal temporel */
    // ndata2=1000*ndata;
    // dates2=span(dates(1),dates(0),ndata2);
    // wd = spl_interp_coefs_new_spline((indgen(0:ndata2-1)/(ndata2-1.))*(ndata-1.), ndata, 3);
    // sig2=wd(spl_spline_samples_to_coefficients(sig,3));
    // avgsig=avg(sig2);
    // sig2=sig2-avgsig;
    // sig=sig-avgsig;
    
    // /* recherche des minimums et maximums locaux */
    // mins=[];
    // maxs=[];
    // idmins=[];
    // idmaxs=[];
    
    // for(k=2;k<=ndata2-1;++k) {
    //     if (sig2(k-1)>sig2(k) & sig2(k+1)>=sig2(k)) {
    //         grow, mins, sig2(k);
    //         grow, idmins, k;
    //     }
    //     else if (sig2(k-1)<sig2(k) & sig2(k+1)<=sig2(k)) {
    //         grow, maxs, sig2(k);
    //         grow, idmaxs, k;
    //     }
    // }
    
    // // window, 1; fma; plg, sig, dates; // plg, sig2, dates2, color="red";
    // // window, 1; plg, zeros, dates2(indices), type="none", marks=1, marker='\4', color="blue";
    
    // /* P�riode moyenne */
    // tcycl = 0.5*(avg(abs((dates2(idmins))(dif)))+avg(abs((dates2(idmaxs))(dif))));
    // half_tcycl = 0.5*tcycl;
    // dates3=dates2(idmins);
    // dates4=dates2(idmaxs);
    // grow, dates4, dates3;
    // dates4=dates4(sort(dates4));

    // /* tronquer les p�riodes incompl�tes */
    // idtronc = where(dates>=dates3(1) & dates<dates3(0));
    // if(is_array(idtronc)) {
    //     dates = dates(idtronc);
    //     theta = theta(idtronc);
    //     data = data(..,idtronc);
    //     u_off = u_off(idtronc);
    //     v_off = v_off(idtronc);
    //     ndata = numberof(theta);
    //     sig3=sig(idtronc);
    // }
    
    // /* Normalisation des dates par les demi-p�riodes moyennes */
    // nextr = numberof(dates4);
    // newdates = array(double,ndata);
    // t_old=dates4(1);
    // for (u=1; u<=nextr-1; ++u) {
    //     Ttemp = dates4(u+1)-dates4(u);
    //     idTtemp = where(dates>=dates4(u) & dates<dates4(u+1));
    //     if(is_array(idTtemp)) {
    //         newdates(idTtemp)=(dates(idTtemp)-t_old)*(half_tcycl/Ttemp)+(u-1)*half_tcycl;
    //     }
    //     t_old = dates4(u+1);
    // }
    
    // // window, 1; plg, sig3, dates3(1)+newdates, color="blue";
    
    // /* temporal spline degree */
    // t_deg = 1;

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

    dweights = array(1.0,nu,nv,ndata);
    id_invalid = where(data < 0.0);
    if (is_array(id_invalid))
        dweights(id_invalid) = 0.0;
    dweights(,1:5,)=0.0;
    dweights(,nv-4:nv,)=0.0;
    dweights(1:5,..)=0.0;
    dweights(nu-4:nu,..)=0.0;

    Htomo_init=yhd_restore("./preliminary_results/CLBfantome_SD_95x95x72x48_voxel_4mm.rec");
    x_init=h_get(Htomo_init,"XR_mus1_mut1_eps1-2").x;

    nx_init=ny_init=95; nz_init=48;
    
    wop = spl_spline_interpolator((indgen(0:nx-1)/(nx-1.))*(nx_init-1.), nx_init, 3, (indgen(0:ny-1)/(ny-1.))*(ny_init-1.), ny_init, 3, (indgen(0:nz-1)/(nz-1.))*(nz_init-1.)., nz_init, 3, indgen(0:nt-1), nt, 0);
    x_init = wop(spl_spline_samples_to_coefficients(x_init,3,3,3,0));

    h_set_copy, Htomo.X, "voxels", x_init;
    
    h_set, Htomo.Cops, dweights=dweights;
    trtt_4D_save, Htomo, Htomo_name, overwrite=1;
    Htomo = trtt_4D_load(Htomo_name);

    /* Create results hash table */
    Htomo_reconst = h_new();
    yhd_save, Htomo_reconst_name, Htomo_reconst, overwrite=1;
}

local Cops; eq_nocopy, Cops, Htomo.Cops;

/* GO RECONSTRUCTION ! */
regulTV = h_new(weight=[mu_s, mu_s, mu_s, mu_t], threshold = eps, options = RGL_TOTVAR_ISOTROPIC);
// regulTV = h_new(weight=[mu_s, mu_t], threshold = [eps1,eps2], flag_separable=1n);

x_iter = Htomo_reconst.x_iter;
if (is_void(x_iter)) x_iter=array(double,Htomo.X.nx,Htomo.X.ny,Htomo.X.nz,Htomo.X.nt);

for (i=1; i<=20; ++i) {
    XR = trtt_4D_optim_simu_launcher(Cops, trtt4D_cost_quadratic_mpy_opky, x=x_iter, mem=5, dweights=dweights, xmin=0.0, regulTV=regulTV, maxiter=10, verbose=1n);
    
    h_set_copy, Htomo_reconst, XRname, XR;
    x_iter = XR.x;
    h_set_copy, Htomo_reconst, "x_iter", x_iter;
    yhd_save, Htomo_reconst_name, Htomo_reconst, overwrite=1;    
}

trtt_4D_save, Htomo, Htomo_name, overwrite=1;

mp_exec, "if (!mp_rank) quit;";
