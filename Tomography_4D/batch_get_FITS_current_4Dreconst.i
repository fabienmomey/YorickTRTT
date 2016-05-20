/*** PETIT SCRIPT POUR ENREGISTRER LA RECONST COURANTE EN FITS (POUR VISU SUR ICY PAR EXEMPLE) ***/
include, "Ytrtt4D.i", 1;
include, "fits.i", 1;
datadir="preliminary_results_from_kubilai/";
Htomo_reconst_name="./preliminary_results_from_kubilai/CLBpatient_SD_120x120x72x13_voxel_8mm_detector_8mm_new.rec";
Htomo_reconst = yhd_restore(Htomo_reconst_name);
// x_iter = Htomo_reconst.x_iter;
x_iter = h_get(Htomo_reconst,"XRglob_mus1e3_mut1e5_eps1e-6").x;
nt=13;
for (i=1;i<=nt;++i){
    fits_write, swrite(format="%sx_detec8mm_1e3_1e5_%02d.fits",datadir,i), x_iter(:,:,:,i), overwrite=1n;
}

