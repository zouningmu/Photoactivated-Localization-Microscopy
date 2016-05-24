function [PC_sgmErflt]=iQ_sgmErflt(DIC,mask_Red_grid,XYZ_count,CAM_set,EM_wave_FP,IM_info,FIL)
%% THIS FUNCTION FILTER THE XYZ_COUNT DATA WITH VARIOUS FILTERS 
% INPUTS:
% DIC: DIC image;
% mask_Red_grid: two column data, xy coordinates inside each cell mask
% XYZ_count: 16 column matrix,fitting results of all spots.
% CAM_set=[CCD_gain, CCD_sens, CCD_QE]; CCD_gain=300;   CCD_sens=67.12; CCD_QE=0.97; 
% EM_wave_FP: emission wavelength of dye
% IM_info=[Pixel_x, Pixel_y]; Pixel_x= 135.8; Pixel_y=176.6; 
% FIL=[ 
%   FIL_sgmx    FIL_sgmy    FIL_ellip   FIL_errx    FIL_erry    FIL_int
%   MIN_sgmx    MIN_sgmy    MIN_ellip   MIN_errx    MIN_erry    MIN_int
%   MAX_sgmx    MAX_sgmy    MAX_ellip   MAX_errx    MAX_erry    MAX_int ]
%   sgmx,sgmy,errx,erry are in nm
%   ellip is unitless, int: counts
% OUTPUTS:
% PC_sgmErflt: cell data, each cell represents one E.coli and contain 16 column matrix after filtering

%% TOP of the routine
FIL_sgmx=FIL(1,1);    FIL_sgmy=FIL(1,2);    FIL_ellip=FIL(1,3);
FIL_errx=FIL(1,4);    FIL_erry=FIL(1,5);    FIL_int=FIL(1,6);
%% physical constants
h = 6.626e-34;    % Planck Constant (m^2 kg / s)
c = 299792458;    % Speed of light (m / s)
ev = 1.602e-19;    % Energy of 1eV (kg m^2 / s^2)
E_ph= h*(c/(EM_wave_FP*1e-9));

%% image information
Pixel_x=IM_info(1); Pixel_y=IM_info(2);

%% find out spots inside the cell mask
imarea=zeros(size(DIC));
imarea(sub2ind(size(imarea),mask_Red_grid(:,2),mask_Red_grid(:,1)))=1;
cc = bwconncomp(imarea, 4);
mask_grid=cell(length(cc.PixelIdxList),1);
for ii=1:length(cc.PixelIdxList);
    grid=false(size(imarea));
    grid(cc.PixelIdxList{ii})=true;
    [row_cc,col_cc] = find(grid>0);
    mask_grid{ii,1}=[col_cc,row_cc];
end
XYZ_count_its=[round(XYZ_count(:,1)), round(XYZ_count(:,2))];
PC=cell(length(mask_grid),1);
for m=1:length(mask_grid)  % loop through all cells
    ptn_member=ismember(XYZ_count_its,mask_grid{m,1},'rows');
    XYZ_coor=XYZ_count(ptn_member,:);
    PC{m,1}=XYZ_coor;                           % protein coordinates inside each cell     
end
%% calculate the Errx and Erry for each spot
PC_sgmErflt=cell(length(PC),1);
for i=1:length(PC);
    if isempty(PC{i,1}); continue;
    else
        PC_mat=PC{i,1};
        cts_n=PC_mat(:,15);
        cts_b=PC_mat(:,16);
        g=CAM_set(1); S=CAM_set(2); QE=CAM_set(3);
        N_e=(cts_n/g)*(S/QE)*3.65;
        N_p=N_e*ev/E_ph;
        b_e=(cts_b/g)*(S/QE)*3.65;
        b_p=b_e*ev/E_ph;

        sgm_x=PC_mat(:,8);   % PC_mat(:,8) gives the sigmax in nm
        Err_x= ((sgm_x.^2)./N_p+...
                Pixel_x^2./(12*N_p)+...
                (8*pi*(sgm_x.^4).*b_p.^2)./(Pixel_x*N_p).^2).^0.5;

        sgm_y=PC_mat(:,9);  % PC_mat(:,9) gives the sigmay in nm
        Err_y= ((sgm_y.^2)./N_p+...
                Pixel_y^2./(12*N_p)+...
                (8*pi*(sgm_y.^4).*b_p.^2)./(Pixel_y*N_p).^2).^0.5;

   % PC_new=[x(px) y(px) x(nm) y(nm) z(nm) 
   %        sigma_x(px) sigma_y(px) sigma_x(nm) sigma_y(nm) ellip PSFW_H(px) 
   %        frame# amp(cts) background(cts) Int(cts) std_b(cts)  Err_x(nm) Err_y(nm)] 
 
PC_new=[PC_mat(:,1:16) Err_x Err_y]; 

%% filter the data
% sigmax filter
if FIL_sgmx == 1; MIN_sgmx=FIL(2,1); MAX_sgmx=FIL(3,1);
    sgmx_raw=PC_new;
    sgmxf= sgmx_raw(:,8)>MIN_sgmx & sgmx_raw(:,8)<MAX_sgmx;
    sgmx_new=sgmx_raw(sgmxf,:);
    flt_new=sgmx_new;
end
% sigmay filter
if FIL_sgmy == 1; MIN_sgmy=FIL(2,2); MAX_sgmy=FIL(3,2);
    sgmy_raw=flt_new;
    sgmyf= sgmy_raw(:,9)>MIN_sgmy & sgmy_raw(:,9)<MAX_sgmy;
    sgmy_new=sgmy_raw(sgmyf,:);
    flt_new=sgmy_new;
end
% Elli filter
if FIL_ellip == 1; MIN_ellip=FIL(2,3); MAX_ellip=FIL(3,3);
    ell_raw=flt_new;
    ell_f=ell_raw(:,10)>MIN_ellip &ell_raw(:,10)<MAX_ellip ;
    ell_new=ell_raw(ell_f,:);
    flt_new=ell_new;
end
% Errx filter
if FIL_errx == 1; MIN_errx=FIL(2,4); MAX_errx=FIL(3,4);
    errx_raw=flt_new;
    errxf= errx_raw(:,17)>MIN_errx & errx_raw(:,17)<MAX_errx;
    errx_new=errx_raw(errxf,:);
    flt_new=errx_new;
end
% Erry filter
if FIL_erry == 1; MIN_erry=FIL(2,5); MAX_erry=FIL(3,5);
    erry_raw=flt_new;
    erryf= erry_raw(:,18)>MIN_erry & erry_raw(:,18)<MAX_erry;
    erry_new=erry_raw(erryf,:);
    flt_new=erry_new;
end
% Int filter
if FIL_int == 1; MIN_int=FIL(2,6); MAX_int=FIL(3,6);
    int_raw=flt_new;
    int_f=int_raw(:,15)>MIN_int &int_raw(:,15)<MAX_int;
    int_new=int_raw(int_f,:);
    flt_new=int_new;
end

PC_new_flt=flt_new;

    end
    PC_sgmErflt{i,1}=PC_new_flt;
end


