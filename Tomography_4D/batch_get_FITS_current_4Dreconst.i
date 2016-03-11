/*** PETIT SCRIPT POUR ENREGISTRER LA RECONST COURANTE EN FITS (POUR VISU SUR ICY PAR EXEMPLE) ***/
include, "Ytrtt4D.i", 1;
include, "fits.i", 1;
datadir="preliminary_results_from_kubilai/";
Htomo_reconst_name="./preliminary_results_from_kubilai/CLBpatient_SD_120x120x72x13_voxel_4mm.rec";
Htomo_reconst = yhd_restore(Htomo_reconst_name);
x_iter = Htomo_reconst.x_iter;
nt=13;
for (i=1;i<=nt;++i){
    fits_write, swrite(format="%sx_iter%02d.fits",datadir,i), x_iter(:,:,:,i), overwrite=1n;
}
