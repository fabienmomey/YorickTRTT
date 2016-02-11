%************************************
%Lecture et affichage des projections
%************************************
clear all;
close all;
clc;

fic='Pj90_shepp.data';
directory='';

nc=1024;
nl=1024;
Nproj=90;

type=1; %1=raw, 2=dcm, 3=bw, 4=tif
switch type
    case 1
        fid = fopen([directory,fic],'r'); %noo
        count=128; %noo
        H = fread(fid, count,'float');
    case 2
        info = dicominfo([directory,fic])
        fid = fopen([directory,fic],'r'); %dcm
        count=3456; %dicom
        H = fread(fid, count,'float32');
    case 3
        fid = fopen([directory,fic],'r'); %bw
        count=128;
        H = fread(fid, count,'float32');
    case 4
end

count=nl*nc;

for i=1:1:Nproj
    fprintf('\nProcessing projection %d',i);

    switch type
        case 1
            I = fread(fid, count,'float'); %raw
        case 2
            I = fread(fid, count,'uint16'); %dcm
        case 3
            I = fread(fid, count,'uint8'); %bw
        case 4
            I=imread([directory,fic],i); %tif
    end

    I = reshape(I,nl,nc);
    figure(1), imagesc(I),colormap gray,title(int2str(i));
    pause(0.1);
end

if type~=4
    fclose(fid);
end


