%% This iQPALM program is to quantitatively count the number of proteins inside the E.coli cell.

%% General Parameter inputs:
% image information
Pixel_x= 135.4; Pixel_y=176.1; 
IM_info=[Pixel_x, Pixel_y];
% Camera setting
CCD_gain=150;   CCD_sens=67.12;  CCD_QE=0.97; 
CAM_set=[CCD_gain, CCD_sens, CCD_QE];
% fluorescence protein emission wavelength
EM_wave_FP=584;

%% Data inputs
dir_DIC='path to directory contains DIC images'; 
cd(dir_DIC); D_DIC=dir;
dir_Red='path to directory contains ensemble red images'; 
cd(dir_Red); D_Red=dir;
dir_Green='path to directory contains ensemble green images'; 
cd(dir_Green); D_Green=dir;

DIC=imread(D_DIC(3).name);  
Red=imread(D_Red(3).name);
Green=imread(D_Green(3).name);

img_zcal='path to directory contains zcal images';
img_dc='path to directory contains red mEos2 images';

%% iQPALM parameter inputs
% for cell boundary and spot finding
dil_sz=3;     threshold=4;    mask_sz=13;   

% repeated spot criteria, Delta=[max_x,max_y,max_z];
Delta=[100,100,200];

% spot filter, use to remove bad spots
% FIL=[ 
%   FIL_sgmx    FIL_sgmy    FIL_ellip   FIL_errx    FIL_erry    FIL_int
%   MIN_sgmx    MIN_sgmy    MIN_ellip   MIN_errx    MIN_erry    MIN_int
%   MAX_sgmx    MAX_sgmy    MAX_ellip   MAX_errx    MAX_erry    MAX_int ]
FIL=[ 1,    1,    0,      1,      1,      0;...
      120,   120,   0,      0,      0,      100;...
      500,  500,  0.15,   1000,   1000,   1000000];

% for reconstructed image
sigma_RI=1 ;     mask_size_RI=11;  precision_RI=10;              
RI_para=[mask_size_RI, sigma_RI,precision_RI];

%% iQPALM main routine
% calibrate the pixel size
[Pixel_x,Pixel_y] = iQ_pixcal(x_cal,y_cal,12.5);
% get z calibration curve
[Z_pos,Z_cal_curve,XYPSF,fitResult,cf_,goodness] = iQ_Zcal(img_zcal,z_step,mask_sz,IM_info);
% find the cell boundary
[DIC_bd,mask_DIC_grid,DIC_coor]=iQ_cellfnd(DIC,dil_sz,savefig);
% correct for the Red image shift
[Red_DIC, mask_Red_grid,Rshf,DCA,RCA ] = iQ_imgshift( DIC,Red,DIC_coor,mask_sz,IM_info,refcoor,Z_cal_curve);
% get drift correction curve
[DFT_pos,DFT_dr,DFT_dr_ave ] = iQ_driftcorr(img_dc,threshold,mask_sz,IM_info,refcoor,Z_cal_curve );
% calculate localization precision using position marker localization
[XYWODC,XYWDC,XZWODC,XZWDC,YZWODC,YZWDC] = iQ_sgmloc( DFT_pos,DFT_dr_ave );
% quatatatively counting the number of proteins inside the cell
[XYZ_count,countprotein,XYZ_repeat,repeatprotein,countall]=iQ_ProtnCount(img_dc,threshold,mask_sz,Delta,mask_Red_grid,DFT_dr_ave,IM_info,Z_cal_curve);
% remove all the spots outside the boundary and allocate spots within the
% same cell into the same group, filter spots with all different kinds of filters and got the correct
% protein number
[PC_sgmErflt]=iQ_sgmErflt(DIC,mask_Red_grid,XYZ_count,CAM_set,EM_wave_FP,IM_info,FIL);
% rotate the cell so that the long axis of the cell lie down along x axis
[celltheta] = iQ_cellrot( DIC, mask_Red_grid,dil_sz,IM_info,PC_sgmErflt );
% profile each cells
[NC,RI,PF_yz,PF_xz,PF_xy,PD]=iQ_profile(PC_sgmErflt,celltheta,Rshf,IM_info,Ausz,RCA,RI_para,bina,binb,binc,savefig);
