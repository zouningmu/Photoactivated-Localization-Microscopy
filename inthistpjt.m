function [ XPJT,YPJT,xx,yy,XY ] = inthistpjt( rimg_x,rimg_y,roi_xranpx,roi_yranpx,binb,binc,px_sz )
%This function generate the 2D intensity histogram and its x and projection
% INPUTS: 
% rimg_x: x data
% rimg_y: y data
% roi_xrange: range of roi in x
% roi_yrange: ronge of roi in y
% binb: number of bin in x
% binc: number of bin in y
% px_sz: pixel size
% OUTPUTS:
% XPJT: two column data, projection of 2D histogram on x
% YPJT: two column data, projection of 2D histogram on y
% xx,yy: mesh grid for XY
% XY: 2D histogram

%% TOP of the routine

xd = rimg_x; yd = rimg_y;
xi = linspace(1,roi_xranpx+1,binb); 
yi = linspace(1,roi_yranpx+1,binc);
xr = interp1(xi,0.5:numel(xi)-0.5,xd,'nearest'); 
yr = interp1(yi,0.5:numel(yi)-0.5,yd,'nearest');
XY = accumarray([yr xr]+0.5, 1, [binb binc]);    
[xx, yy] = meshgrid...
    (roi_xranpx/binb:roi_xranpx/binb:roi_xranpx,...
     roi_yranpx/binc:roi_yranpx/binc:roi_yranpx);
xx = xx*px_sz; yy = yy*px_sz;

XY_xpjt_ax = xx(1,:)+0.5*xx(1,1);  XY_xpjt = sum(XY); 
XY_ypjt_ax = yy(:,1)+0.5*yy(1,1);  XY_ypjt = sum(XY,2);


XPJT=[XY_xpjt_ax',XY_xpjt'];
YPJT=[XY_ypjt_ax',XY_ypjt'];
end

