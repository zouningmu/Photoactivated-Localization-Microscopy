 function [ gf_xypsfw_h ] = iQ_gaussianfit( pick,im_bk,mask_sz,n,IM_info,Z_cal_curve )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
Pixel_x=IM_info(1); Pixel_y=IM_info(2);  corr_cylin=IM_info(1)/IM_info(2);

gf_xypsfw_h=cell(1,length(pick)/2);
for j=1:length(pick)/2; 
   [xi,yi] = meshgrid(pick(1,2*j-1)-floor(mask_sz/2):pick(1,2*j-1)+floor(mask_sz/2)...
       ,pick(1,2*j)-floor(mask_sz/2):pick(1,2*j)+floor(mask_sz/2));
   zi=im_bk(pick(1,2*j)-floor(mask_sz/2):pick(1,2*j)+floor(mask_sz/2)...
       ,pick(1,2*j-1)-floor(mask_sz/2):pick(1,2*j-1)+floor(mask_sz/2));
   results = iQ_autoGaussianSurfML(xi,yi,zi);
   amp=results.a; background=results.b; sigma_x=results.sigmax; sigma_y=results.sigmay/corr_cylin;
   mu_x=results.x0; mu_y=results.y0/corr_cylin; Int=sum(sum(results.G-results.b)); std_b=mean(std(zi-results.G));
   ellip=2*abs((sigma_x-sigma_y))/(sigma_x+sigma_y);
   PSFW_H=sigma_x-sigma_y;     % calculate the width-height vaue from PSF fitting
   mu_z=0;
   
   if nargin==6;
   x=Z_cal_curve(:,1);                                     
   y=Z_cal_curve(:,2);
% calulate z from the width-height of psf
   mu_z=interp1(y,x,PSFW_H);
   end
   
   gf_xypsfw_h{1,j}=[mu_x mu_y mu_x*Pixel_x mu_y*Pixel_x mu_z...
                     sigma_x sigma_y sigma_x*Pixel_x sigma_y*Pixel_x ellip...
                     PSFW_H n amp background Int std_b ];

end

end