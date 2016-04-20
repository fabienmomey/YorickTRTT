include, "Ytrtt2D_dyn.i", 1;

datadir="dataCLB_new_14_01_2016/";
datafile="results_dataset10_380x380x22_pixel1mm";
imagefile="Reconstruction_Grid_MultipleHyperParameters.fits";
paramfile="Reconstruction_Grid_MultipleHyperParameters.txt";
f=create(datadir+paramfile);

/* Dimensions objet 2D+t */
nx = 380;
ny = 380;
nt = 22;

/* Recharger les données de reconstruction */
Htomo=trtt_2D_load(datadir+datafile);
// KEYS=["R_eps1_1e-6_eps2_1e-6",
//       "R_eps1_1e-1_eps2_1e-6",
//       "R_eps1-1_eps2-1e-6",
//       "Reps1e-6",
//       "Reps1e-2",
//       "Reps1e-1",
//       "Reps1e0"];
KEYS=["R_eps1_1e-6_eps2_1e-6",
      "Reps1e-6"];
// KEYS=["R_eps1_1e-6_eps2_1e-6"];

/* h_keys(Htomo) */
// Nbtot=0;
// for (k=1; k<=numberof(KEYS); ++k) {
//     SUBKEYS=h_keys(h_get(Htomo,KEYS(k)));
//     write, format="%s : \t Nb=%d\n", KEYS(k), numberof(SUBKEYS);
//     Nbtot+=numberof(SUBKEYS);
// }
// write, format=" \t Nb total = %d\n", Nbtot;

/* Définir la taille NxN de la grille des ROIs */
Ngrid=15;
/* Définir la ROI à extraire (position et dimensions) */
xroi=214;
yroi=174;
Nroi=64;
/* Taille de l'image globale de la grille des ROIs */
Nimage=(Nroi+2)*Ngrid;
Image=array(1.0,Nimage,Nimage,nt);


cpt=0;
for (k=1; k<=numberof(KEYS); ++k) {
    KEY=h_get(Htomo,KEYS(k));
    KEY_keys=h_keys(KEY);
    NKEY=numberof(KEY_keys);
    for (i=1;i<=NKEY;++i) {
        cpt++;
        
        /* Récupérer les inforamtions de la reconstruction courante */
        XR=h_get(KEY,KEY_keys(i));
        weight=XR.regulTV.weight;
        threshold=XR.regulTV.threshold;
        flag_separable=XR.regulTV.flag_separable;
        /* Extraire ROI */
        roi=XR.x(xroi:(xroi+Nroi-1),yroi:(yroi+Nroi-1),);
  
        /* Positionner la ROI dans l'image globale */
        m=long((cpt-1)%Ngrid+1);
        l=long(floor((cpt-1)/Ngrid)+1);

        Image((m-1)*(Nroi+2)+1:(m-1)*(Nroi+2)+Nroi,(l-1)*(Nroi+2)+1:(l-1)*(Nroi+2)+Nroi,)=roi;

        /* Écrire les infos des hyperparamètres dans le fichier texte */
        if (!is_void(flag_separable))
            write, f, format="N°%d \t ; \t Case [%d,%d] \t ; \t REGUL SEPA \t ; \t threshold=[%.2e,%.2e] \t ; \t weight=[%.2e,%.2e]\n", cpt, m, l, threshold(1), threshold(2), weight(1), weight(2);
        else
            write, f, format="N°%d \t ; \t Case [%d,%d] \t ; \t REGUL GLOB \t ; \t threshold=%.2e \t ; \t weight=[%.2e,%.2e]\n", cpt, m, l, threshold, weight(1), weight(3);
    }
}
/* Enregistrer l'image globale */
fits_write, datadir+imagefile, Image, overwrite=1n;
n=close(f);

                                       
