mp_include, "Ytrtt4D.i";
// include, "Ytrtt4D.i", 1;

Htomo_name="./preliminary_results/CLBpatient_SD_120x120x72x13_voxel_4mm_detector_4mm";
Htomo_reconst_name="./preliminary_results/CLBpatient_SD_120x120x72x13_voxel_4mm_detector_4mm.rec";

mu_s = 1.e3;
mu_t = 1.e3;
eps = 1.e-6; // eps1=1.0; eps2=1.e-6;
XRname = "XRglob_mus1e3_mut1e3_eps1e-6"; // XRname = "XRsepa_mus1e0_mut1e0_eps1_1e0_eps2_1e-6";

if (!is_void(open(Htomo_name, "rb", 1))) {
    write, format="---------- %s reloaded ----------\n", Htomo_name;
    Htomo = trtt_4D_load(Htomo_name);
    Htomo_reconst = yhd_restore(Htomo_reconst_name);
    local dweights; eq_nocopy, dweights, Htomo.Cops.dweights;
    s_deg = Htomo.X.s_deg;
    t_deg = Htomo.X.t_deg;
} else {
    write, format="---------- %s created ----------\n", Htomo_name;
    /* PATIENT */
    data_dir= "/home/momey/Data/data_CLB_patient_02-07-2012/";
    // data_dir= "/home/momey/Recherche_Tomographie/Data/data_CLB_patient_02-07-2012/img_1.3.46.423632.141000.1169042526.68/";
    /* FANTÔME MÉCANIQUE */
    // data_dir= "/home/momey/Data/data_CLB_09-11-2011/mvt3d/";
    
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
    s_deg = 3;
    size_footprint = 10000;
    
    /* TOMDATA */
    nv = DATA.nv/4;
    nu = DATA.nu/4;
    u_off = DATA.uoff(5:-1);
    v_off = DATA.voff(5:-1); //FIXME: in mm
    u_scl = 4.0*DATA.u_scl; //FIXME: in mm
    v_scl = 4.0*DATA.v_scl; //FIXME: in mm
    t_index = 0.0;

    /* UNDERSAMPLING OF PROJECTIONS BY 4 */
    data=yeti_convolve(data,kernel=[0.,0.25,0.25,0.25,0.25],which=[1,2])(::4,::4,:);
    
    /* PROJECTOR */
    mode = TRTT_CONE_BEAM;
    type = TRTT_SPLINE_DRIVEN;
    
    /* Noise features */                  
    photon_flux = 1.e6; //FIXME: SNR² ~ photon_flux
    read_noise = 5.;
    var_noise_uniform = 1.;
    
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
    cpt=0;
    for (u=1; u<=nextr-1; ++u) {
        Ttemp = dates4(u+1)-dates4(u);
        idTtemp = where(dates>=dates4(u) & dates<dates4(u+1));
        if(is_array(idTtemp)) {
            newdates(idTtemp)=(dates(idTtemp)-t_old)*(half_tcycl/Ttemp)+cpt*half_tcycl;
            cpt++;
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

    dweights = array(1.0,nu,nv,ndata);
    id_invalid = where(data < 0.0);
    if (is_array(id_invalid))
        dweights(id_invalid) = 0.0;
    dweights(,1:5,)=0.0;
    dweights(,nv-4:nv,)=0.0;
    dweights(1:5,..)=0.0;
    dweights(nu-4:nu,..)=0.0;
    
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

for (i=1; i<=5; ++i) {
    XR = trtt_4D_optim_simu_launcher(Cops, trtt4D_cost_quadratic_mpy_opky, x=x_iter, mem=5, dweights=dweights, xmin=0.0, regulTV=regulTV, maxiter=100, verbose=1n);
    
    h_set_copy, Htomo_reconst, XRname, XR;
    x_iter = spl_spline_samples_to_coefficients(XR.x,s_deg,s_deg,s_deg,t_deg);
    h_set_copy, Htomo_reconst, "x_iter", x_iter;
    yhd_save, Htomo_reconst_name, Htomo_reconst, overwrite=1n;    
}

trtt_4D_save, Htomo, Htomo_name, overwrite=1;

mp_exec, "if (!mp_rank) quit;";

/*** TEST EXTENT ***/
// nx = 120;
// ny = 120;
// nz = 72;
// nt = 13;
// for (i=1;i<=100;++i) {
//     x=random(nx,ny,nz,nt);
//     y=random(nx+2,ny+2,nz+2,nt+1);
//     xCTy=sum((x*trtt_4D_extended_object_operator(y,1))(*));
//     Cxy=sum((trtt_4D_extended_object_operator(x)*y)(*));
//     write, format="xCTy=%f \t Cxy=%f \n", xCTy, Cxy;
// }

