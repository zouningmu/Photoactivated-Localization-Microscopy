function [XYZ_count,countprotein,XYZ_repeat,repeatprotein,countall]=iQ_ProtnCount(folderName,threshold,mask_sz,Delta,mask_Red_grid,DFT_dr_ave,IM_info,Z_cal_curve)
%% THE FUNCTION QUANTITATIVELY COUNT THE SPOTS INSIDE THE BOUNDARY OF THE CELL
% INPUTS:
% folderName: name of the foler, has to be a string
% threshold: minimum number of sigma needed for a peak picking, suggested value is 4. 
% mask_sz: size of mask that cover each spots, recommend size:11
% Delta: values to determine whether x1 y1 z1 = x2 y2 z2 or not. The format should be like [x y z], Note: x,y are in nm and z is in nm.
% mask_Red_grid: grid coordinates of Red mask, use to determine whether spots are inside the cell or not.
% DFT_dr_ave: drift correction data, 2 column data
% IM_info: contains pixel size information  and correction factor for cylindrical lens;
% Z_cal_curve: fitting curve of Z calibration

% OUTPUTS:
% XYZ_count:  [mu_x(px) mu_y(px) mu_x(nm) mu_y(nm) mu_z(nm) sigma_x(px) sigma_y(px) sigma_x(nm) sigma_y(nm) ellip PSFW_H(px) n amp background Int std_b] of counted spots
% countprotein: total number of counted spots
% XYZ_repeat: [mu_x(px) mu_y(px) mu_x(nm) mu_y(nm) mu_z(nm) sigma_x(px) sigma_y(px) sigma_x(nm) sigma_y(nm) ellip PSFW_H(px) n amp background Int std_b] of repeated spots
% repeatprotein: total number of repeated spots
% countall: total number of spots
%% TOP of the routine
tic
cd(folderName); D=dir;
num_frame=length(D)-2;

% get the meshgrid size from mask_sz
sz=(mask_sz-1)/2;
% get the values to judge the spots are the same one or not between two connected frames.
del_x=Delta(1); del_y=Delta(2);                         % 2D
if nargin == 8; del_z=Delta(3);end;                     % 3D
% get the pixel and corr_cylin information
Pixel_x=IM_info(1); Pixel_y=IM_info(2);  corr_cylin=IM_info(1)/IM_info(2);
% make Z calibration curve ready to use
if nargin == 8; x=Z_cal_curve(:,1); y=Z_cal_curve(:,2); end;  % 3D
%%  OPEN 1st IMAGE AND GENERATE X1Y1Z1
     n=1;
   while 1;
     im=double(imread(D(n+2).name));                        % load first frame   
     [im_bk,peak_fnd,~,~]=iQ_pkfnd(im,threshold,mask_sz);   % find out the local maximum position in the current frame
     peak=intersect(peak_fnd,mask_Red_grid,'rows');         % pick up only those spots are inside the cell.

     if isempty(peak)
         n=n+1; disp('nothing inside the cell in current frame')
         continue
     else
     XYZ_current=cell(length(peak),1);
  
     % 2D Gaussian fit all spots inside the cell
     [row_peak ~]=size(peak);
        for j=1:row_peak;
            [xi,yi]=meshgrid(peak(j,1)-sz:peak(j,1)+sz,peak(j,2)-sz:peak(j,2)+sz);
            zi=im_bk(peak(j,2)-sz:peak(j,2)+sz,peak(j,1)-sz:peak(j,1)+sz);
            results = iQ_autoGaussianSurfML(xi,yi,zi);
            amp=results.a; background=results.b; sigma_x=results.sigmax; sigma_y=results.sigmay/corr_cylin;
            mu_x=results.x0; mu_y=results.y0; Int=sum(sum(results.G-results.b)); std_b=mean(std(zi-results.G));
            ellip=2*abs((sigma_x-sigma_y))/(sigma_x+sigma_y);
            PSFW_H=sigma_x-sigma_y;     % calculate the width-height vaue from PSF fitting
           
            % if you are working w/ 2D data, this will be excuted.     
            if nargin < 8;
            mu_z=0;
            end
            
            % if you are working w/ 3D data, this will be excuted.
            if nargin == 8;
            % calulate z from the width-height of psf
            mu_z=interp1(y,x,PSFW_H);
            end;
            XYZ_current{j,1}=[mu_x mu_y mu_x*Pixel_x mu_x*Pixel_y mu_z...
                     sigma_x sigma_y sigma_x*Pixel_x sigma_y*Pixel_x ellip PSFW_H ...
                     n amp background Int std_b ];
        end
     break;  
     end;
   end
countall=length(peak);
XYZ_repeat=cell(num_frame,1);
XYZ_count=cell(num_frame,1);

%%  REPEAT THE SAME PROCESS AS 1st IMAGE THROUGH THR MV
for i=n+3:num_frame+2;
     xyz_current=cell2mat(XYZ_current);
     
     im_nxt=double(imread(D(i).name)+1);                     % load the next frame
     [im_bk_nxt,pk_nxt_fnd,~,~]=iQ_pkfnd(im_nxt,threshold,mask_sz); % find out the local maximum positions in the next frame
     pk_nxt=intersect(pk_nxt_fnd,mask_Red_grid,'rows');     % pick up only those spots are inside the cell.
     
     if isempty(pk_nxt)
          XYZ_nxt={};
          XYZ_nxt{1,1}=[0 0 0 0 0 0 0 0 0 0 0 i-2 0 0 0 0];    
          display(['Frame no. ',int2str(i-2),'   Sptos fitted: ',' nothing inside the cell in next frame'])
         
     else
     countall=countall+length(pk_nxt);
     XYZ_nxt=cell(length(pk_nxt),1);
   
     [row_nxt ~]=size(pk_nxt);
     for k=1:row_nxt;
        [xi_nxt,yi_nxt]=meshgrid(pk_nxt(k,1)-sz:pk_nxt(k,1)+sz,pk_nxt(k,2)-sz:pk_nxt(k,2)+sz);
        zi_nxt=im_bk_nxt(pk_nxt(k,2)-sz:pk_nxt(k,2)+sz,pk_nxt(k,1)-sz:pk_nxt(k,1)+sz);
        results_nxt = iQ_autoGaussianSurfML(xi_nxt,yi_nxt,zi_nxt);
        amp_nxt=results_nxt.a; background_nxt=results_nxt.b; sigma_x_nxt=results_nxt.sigmax; sigma_y_nxt=results_nxt.sigmay/corr_cylin;
        mu_x_nxt=results_nxt.x0; mu_y_nxt=results_nxt.y0; Int_nxt=sum(sum(results_nxt.G-results_nxt.b)); std_b_nxt=mean(std(zi_nxt-results_nxt.G));
        ellip_nxt=2*abs((sigma_x_nxt-sigma_y_nxt))/(sigma_x_nxt+sigma_y_nxt);
        PSFW_H_nxt=sigma_x_nxt-sigma_y_nxt;              % calculate the width-height vaue from PSF fitting
        
        % if you are working w/ 2D data, this will be excuted.
        if nargin < 8;
        mu_z_nxt=0;
        end
        
        % if you are working w/ 3D data, this will be excuted.
        if nargin == 8;
        mu_z_nxt=interp1(y,x,PSFW_H_nxt);                             % calulate mu_z from the width-height of psf   
        end
        
        % drift correction in x,y,z direction     
        mu_x_nxt=mu_x_nxt-DFT_dr_ave(i-2,1);
        mu_y_nxt=mu_y_nxt-DFT_dr_ave(i-2,2);
        mu_z_nxt=mu_z_nxt-DFT_dr_ave(i-2,5);

        XYZ_nxt{k,1}=[mu_x_nxt mu_y_nxt mu_x_nxt*Pixel_x mu_y_nxt*Pixel_y mu_z_nxt...
            sigma_x_nxt sigma_y_nxt sigma_x_nxt*Pixel_x sigma_y_nxt*Pixel_x ellip_nxt PSFW_H_nxt...
            i-2 amp_nxt background_nxt  Int_nxt std_b_nxt];
     end
          display(['Frame no. ',int2str(i-2),'   Sptos fitted: ',int2str(row_nxt), '  spots'])
     end
     
      r2=1;  XYZ_rpt=cell(length(XYZ_current),1);
      [row_XYZ_current ~]=size(cell2mat(XYZ_current));
      for k=1:row_XYZ_current;
         x_current=XYZ_current{k,1}(1,3);
         y_current=XYZ_current{k,1}(1,4);
         z_current=XYZ_current{k,1}(1,5);
         
         [row_XYZ_nxt ~]=size(cell2mat(XYZ_nxt));
         for m=1:row_XYZ_nxt;
             x_nxt=XYZ_nxt{m,1}(1,3);
             y_nxt=XYZ_nxt{m,1}(1,4);
             z_nxt=XYZ_nxt{m,1}(1,5);
             % if you are working w/ 2D data, this will be excuted.
             if nargin < 8,
                if abs(x_current-x_nxt)<del_x && abs(y_current-y_nxt)<del_y;   
                XYZ_rpt{r2,1}=XYZ_current{k,1}(1,:); r2=r2+1;
                end
             end
             % if you are working w/ 3D data, this will be excuted.
             if nargin == 8;
                if abs(x_current-x_nxt)<del_x && abs(y_current-y_nxt)<del_y && abs(z_current - z_nxt) <del_z;  % note: for x, y and z: 60 means 60 nm
                XYZ_rpt{r2,1}=XYZ_current{k,1}(1,:); r2=r2+1;
                end
             end
         end
      end
      
      % summary the coordinates for repeated spots      
      xyz_rpt=cell2mat(XYZ_rpt);
      if XYZ_current{1,1}(1,1)>0;
         XYZ_repeat{i-3,1}=xyz_rpt;
      end
      % get the coordinates for no-repeated spots, which is counted proteins   
      if XYZ_current{1,1}(1,1)>0;
         xyz_count=xyz_current(~ismember(xyz_current,xyz_rpt,'rows'),:);
         XYZ_count{i-3,1}=xyz_count;
      end
      % make coordinates of next frame as current frame for continuous comparison
      XYZ_current=XYZ_nxt;
end

    % calculate the number repeated spots
    XYZ_repeat=cell2mat(XYZ_repeat);
    [repeatprotein ~]=size(XYZ_repeat);
    
    % take care of the last frame    
    XYZ_current=cell2mat(XYZ_current);
   if XYZ_current(1,1)>0;
      XYZ_count{num_frame,1}=XYZ_current;
   end
    % calculate the number of proteins
    XYZ_count=cell2mat(XYZ_count);
    XYZ_count=XYZ_count(~any(isnan(XYZ_count),2),:);
    [countprotein ~]=size(XYZ_count);
toc
%% nested function iQ_pkfnd
end