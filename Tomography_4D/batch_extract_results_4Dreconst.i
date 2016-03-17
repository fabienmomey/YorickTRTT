include, "Ytrtt4D.i", 1;

datadir="preliminary_results/";
datafile="CLBpatient_SD_120x120x72x13_voxel_4mm.rec";
imagefile="Reconstruction_Grid_MultipleHyperParameters.fits";
paramfile="Reconstruction_Grid_MultipleHyperParameters.txt";
f=create(datadir+paramfile);

/* Dimensions objet 2D+t */
nx = 120;
ny = 120;
nz = 72;
nt = 13;

/* Recharger les données de reconstruction */
Htomo=yhd_restore(datadir+datafile);
KEYS=h_keys(Htomo);
// Nbtot=0;
// for (k=1; k<=numberof(KEYS); ++k) {
//     SUBKEYS=h_keys(h_get(Htomo,KEYS(k)));
//     write, format="%s : \t Nb=%d\n", KEYS(k), numberof(SUBKEYS);
//     Nbtot+=numberof(SUBKEYS);
// }
// write, format=" \t Nb total = %d\n", Nbtot;

/* Définir la taille NxN de la grille des ROIs */
Ngrid=2;
/* Définir la ROI à extraire (position et dimensions) */
xroi=7;
yroi=24;
zroi=1;
Nroi=75;
Nroi_z=nz;
/* Taille de l'image globale de la grille des ROIs */
Nimage=(Nroi+2)*Ngrid;
Nimage_z=Nroi_z;
Image=array(1.0,Nimage,Nimage,Nimage_z,nt);

for (k=1; k<=numberof(KEYS); ++k) {
    KEY=h_get(Htomo,KEYS(k));
        
    /* Récupérer les inforamtions de la reconstruction courante */
    XR=h_get(Htomo,KEY);
    weight=XR.regulTV.weight;
    threshold=XR.regulTV.threshold;
    flag_separable=XR.regulTV.flag_separable;
    /* Extraire ROI */
    roi=XR.x(xroi:(xroi+Nroi-1),yroi:(yroi+Nroi-1),zroi:(zroi+Nroi_z-1),);
    
    /* Positionner la ROI dans l'image globale */
    m=k%Ngrid;
    l=long(floor((k-1)/Ngrid)+1);
    
    Image((m-1)*(Nroi+2)+1:(m-1)*(Nroi+2)+Nroi,(l-1)*(Nroi+2)+1:(l-1)*(Nroi+2)+Nroi,1:Nroi_z,)=roi;
    
    /* Écrire les infos des hyperparamètres dans le fichier texte */
    if (!is_void(flag_separable))
        write, f, format="N°%d \t ; \t Case [%d,%d] \t ; \t REGUL SEPA \t ; \t threshold=[%.2e,%.2e] \t ; \t weight=[%.2e,%.2e]\n", k, m, l, threshold(1), threshold(2), weight(1), weight(2);
    else
        write, f, format="N°%d \t ; \t Case [%d,%d] \t ; \t REGUL GLOB \t ; \t threshold=%.2e \t ; \t weight=[%.2e,%.2e]\n", k, m, l, threshold, weight(1), weight(4);
}
/* Enregistrer l'image globale */
fits_write, datadir+imagefile, Image, overwrite=1n;
n=close(f);
