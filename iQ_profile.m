function [NC,RI,PF_yz,PF_xz,PF_xy,PD]=iQ_profile(PC_sgmErflt,celltheta,Rshf,IM_info,Ausz,RCA,RI_para,bina,binb,binc,savefig)
%% THIS FUNCTION ANALYZE THE PROFILES ALONG X,Y,AND Z AXIS
% INPUT: 
% PC_sgmErflt: Protein count data after filtering
% celltheta: contains angle info to rotate a cell along x axis

% Rshf: 
% Rshf=[Red_xshf,Red_yshf,mean(DIC_Au_x-Red_Au_x)+dx_R,mean(DIC_Au_y-Red_Au_y)+dy_R];
% 4 element array, Rshf(1) and (2) give the rounded image shift x and y used for cell mask shifting
%                  Rshf(3) and (4) give the exact image shift x,y 

% DCA: unit: pixel
% DCA(k,1)=xbar; DCA(k,2)=ybar;
% 2 column matrix, column 1 and 2 give the centroid position of kth reference point 

% RCA: unit: pixel
% RCA=[mu_x,mu_y,sigma_x,sigma_y,PSFW_H, mu_z];
% 5 column matrix, PSF fitting result of reference point, column 1 to
% 5 give ceter x,y, sigmax,sigmay,PSF width-height and correspoding z
% position. Note: here mu_y and sigma_y have been corrected for cylindrical distortion.

% IM_info: px information
% Ausz: diameter of the Au marker.
% RI_para: reconstructed image parameter,[ mask_size, sigma,precision]
% bina: number of slices along axis
% binb: grid size
% binc: grid size
% savefig: save figure file or not, input 1 if you want to save a figure or
%          do not input any thing so that the fig wont be saved.
% OUTPUT:
% NC: normalized cell coordinates
% RI: reconstructed image
% PF: profile data
% PD: protein density

%%
close all;
mask_size=RI_para(1); sigma=RI_para(2); precision=RI_para(3);
px_sz=IM_info(1); py_sz=IM_info(2);
delx=(Rshf(1)-Rshf(3))*px_sz; % in nm
dely=(Rshf(2)-Rshf(4))*py_sz; % in nm
zau=RCA(6);             % in nm 
rau=Ausz/2;             % in nm

if nargin<11;    savefig=0; end;
if savefig==1; 
%     dirname=input('input ppt name: ','s');
%     pptfilname=strcat(dirname,'.ppt');
      pptfilname=('cellprofile_RI.ppt');
      pptfilname1=('cellprofile_RIslince.ppt');
      pptfilname2=('cellprofile_Intslince.ppt');
end;
%% ROTATE THE CELL IN XY PLANE
% find out the rotation angles for every cells using cellmask
  [r_celltheta,~]=size(celltheta);
  RI = cell(r_celltheta,5);
  PF_yz = cell(r_celltheta,bina);
  PF_xz = cell(r_celltheta,bina);
  PF_xy = cell(r_celltheta,bina);
  PD = cell(r_celltheta,4);
  NC = cell(r_celltheta,3);
  for ii=1:r_celltheta
      close all
 
      PC_sgmErflt_mat=PC_sgmErflt{ii};
      [row_sgmErflt,~]=size(PC_sgmErflt_mat);
      if row_sgmErflt < 2;
          continue
      else
          
      theta_cm = cell2mat(celltheta(ii,1));
      rec_para = cell2mat(celltheta(ii,2));
      cell_vol = cell2mat(celltheta(ii,3));
          
      lcx=rec_para(1); lcy=rec_para(2); lenc=rec_para(3);  widc=rec_para(4);
      lencpx=lenc/px_sz; widcpx=widc/px_sz;
      
      xxx=PC_sgmErflt_mat(:,3); yyy=PC_sgmErflt_mat(:,4); zzz=PC_sgmErflt_mat(:,5);  % in nm
      % correct for imgshift with better resolution
      xxx=xxx-delx; yyy=yyy-dely;
      
      ProtnNum = length(xxx);
      ProtnDen = round(ProtnNum/cell_vol*10)/10;
  %------------------------------------------------------------------------    
      PD{ii,1}=rec_para; PD{ii,2}=ProtnNum; PD{ii,3}=cell_vol; PD{ii,4}=ProtnDen;
  %----------------------------------------------------------------------    
      % convert cartesian coordinates to spherical coordinates
      r=(xxx.^2+yyy.^2).^0.5;
      t=atan(yyy./xxx);
      t2=t-theta_cm;
      xxx_new=r.*cos(t2);
      yyy_new=r.*sin(t2);
      
      % move the axis 
      xxx=xxx_new-lcx+widc/2+px_sz;  % in nm
      yyy=yyy_new-lcy+widc/2+px_sz;  % in nm
      zzz=zzz-zau+rau+widc/2+px_sz;  % in nm   
      
      %ROI range
      % in nm
      roi_xmin = px_sz; roi_xmax = lenc+widc+px_sz; roi_xran=roi_xmax-roi_xmin;
      roi_ymin = px_sz; roi_ymax = widc+widc+px_sz; roi_yran=roi_ymax-roi_ymin;
      roi_zmin = px_sz; roi_zmax = widc+widc+px_sz; roi_zran=roi_zmax-roi_zmin;
      % in px
      roi_xminpx = roi_xmin/px_sz; roi_xmaxpx = roi_xmax/px_sz; roi_xranpx = roi_xran/px_sz;
      roi_yminpx = roi_ymin/px_sz; roi_ymaxpx = roi_ymax/px_sz; roi_yranpx = roi_yran/px_sz;
      roi_zminpx = roi_zmin/px_sz; roi_zmaxpx = roi_zmax/px_sz; roi_zranpx = roi_zran/px_sz;
      
      % Remove spots outside the ROI (in nm)
      xyzdata=[xxx,yyy,zzz];
      ix= xyzdata(:,1)>2*roi_xmin & xyzdata(:,1) < roi_xmax;
      xyzdata=xyzdata(ix,:);
      iy= xyzdata(:,2)>2*roi_ymin & xyzdata(:,2) < roi_ymax;
      xyzdata=xyzdata(iy,:);
      iz= xyzdata(:,3)>2*roi_zmin & xyzdata(:,3) < roi_zmax;
      xyzdata=xyzdata(iz,:);
      xxx=xyzdata(:,1); yyy=xyzdata(:,2); zzz=xyzdata(:,3);
      
      % Red data, convert nm into px
      rimg_x=xxx/px_sz;   % in pixel  
      rimg_xstep=roi_xranpx/bina;
     
      rimg_y=yyy/px_sz;   % in pixel
      rimg_ystep=roi_yranpx/bina;
      
      rimg_z=zzz/px_sz;   % in pixel
      rimg_zstep=roi_zranpx/bina; 
      
      % Normalized x,y,z : the coordinates are normalized to a cell shape with aspect ratio of 1:2 for width and length respectively.
      % 0 represents the left and bottom of the cell, 1 and 2 represents the limit of the width and length
      nrimg_x=(rimg_x-(widcpx/2+1))/(lencpx/2);
      nrimg_y=(rimg_y-(widcpx/2+1))/widcpx;
      nrimg_z=(rimg_z-(widcpx/2+1))/widcpx;
%--------------------------------------------------------------------------      
      NC{ii,1}=nrimg_x; NC{ii,2}=nrimg_y; NC{ii,3}=nrimg_z;
%--------------------------------------------------------------------------      
      % here we +1 in the roi_range due to the starting point is already 1 px away from origin.     
      RI_yz=iQ_ReconImage([rimg_y,rimg_z],[roi_yranpx+1,roi_zranpx+1],mask_size,sigma,precision);
      RI_xz=iQ_ReconImage([rimg_x,rimg_z],[roi_xranpx+1,roi_zranpx+1],mask_size,sigma,precision);
      RI_xy=iQ_ReconImage([rimg_x,rimg_y],[roi_xranpx+1,roi_yranpx+1],mask_size,sigma,precision);
      
      % generate intensity histogram for yz xz and xy
      [ YZ_xpjt,~,YZ_yy,YZ_zz,YZ ] = inthistpjt( rimg_y,rimg_z,roi_yranpx,roi_zranpx,binb,binc,px_sz );
      [ XZ_xpjt,~,XZ_xx,XZ_zz,XZ ] = inthistpjt( rimg_x,rimg_z,roi_xranpx,roi_zranpx,binb,binc,px_sz );
      [ XY_xpjt,~,XY_xx,XY_yy,XY ] = inthistpjt( rimg_x,rimg_y,roi_xranpx,roi_yranpx,binb,binc,px_sz );
%------------------------------------------------------------------------     
      RI{ii,1}=RI_xy; RI{ii,2}=RI_xz; RI{ii,3}=RI_yz; 
      RI{ii,4}=[rimg_x,rimg_y,rimg_z];
      RI{ii,5}=[rimg_xstep,rimg_ystep,rimg_zstep,roi_xranpx,roi_yranpx,roi_zranpx];
%------------------------------------------------------------------------     
      currentdic = pwd;
      fig_RI=figure; 
      ah1 = subplot(1,5,1);
        imshow(RI_yz,[]); set(ah1,'ydir','normal'); colormap hot; title('yz int hist');
        rectangle('Position',[(widcpx/2+1)*10,(widcpx/2+1)*10,widcpx*10,widcpx*10],'Curvature',[1,1],'LineWidth',2,'LineStyle','--','EdgeColor','w');
        axis([roi_yminpx*10,roi_ymaxpx*10,roi_zminpx*10,roi_zmaxpx*10]); 
      ah4 = subplot(1,5,[2,3]); 
        imshow(RI_xz,[]); set(ah4,'ydir','normal'); colormap hot; title('xz int hist');
        rectangle('Position',[(widcpx/2+1)*10,(widcpx/2*10+1),lencpx*10,widcpx*10],'Curvature',1,'LineWidth',2,'LineStyle','--','EdgeColor','w');
        axis([roi_xminpx*10,roi_xmaxpx*10,roi_zminpx*10,roi_zmaxpx*10]); 
      ah7 = subplot(1,5,[4,5]); 
        imshow(RI_xy,[]); set(ah7,'ydir','normal'); colormap hot; title('xy int hist');
        rectangle('Position',[(widcpx/2+1)*10,(widcpx/2+1)*10,lencpx*10,widcpx*10],'Curvature',1,'LineWidth',2,'LineStyle','--','EdgeColor','w');
        axis([roi_xminpx*10,roi_xmaxpx*10,roi_yminpx*10,roi_ymaxpx*10]);
      
      if savefig==1;
      saveppt2(pptfilname,'figure',fig_RI,'scale',true,'stretch',false,'editmenu','title',['Cell no. ' num2str(ii),'Protn Den:' num2str(ProtnDen) ' um-3']);
      end
        
      fig_RI2=figure;  
      ah2 = subplot(2,5,1);
        pcolor(YZ_yy,YZ_zz,YZ); cbh=colorbar('East'); colormap jet; set(cbh,'Ycolor','r')
        rectangle('Position',[widc/2+px_sz,widc/2+px_sz,widc,widc],'Curvature',[1,1],'LineWidth',2,'LineStyle','--','EdgeColor','w');
        axis equal; axis([roi_ymin,roi_ymax,roi_zmin,roi_zmax]); 
        set(ah2,'XTickLabel','','YTickLabel',''); 
      ah3 = subplot(2,5,6);
        nYZ_xpjt_x=(YZ_xpjt(:,1)-((widcpx/2+1)*px_sz))/widc;
        bar(ah3,nYZ_xpjt_x,YZ_xpjt(:,2),1,'g');
        nlim_min = (roi_ymin-(widc/2+px_sz))/widc;
        nlim_max = (roi_ymax-(widc/2+px_sz))/widc;
        xlim([nlim_min,nlim_max]); xlabel('norm position'); ylabel('Counts');
        hold on;
        ax3t =  axes('Position',get(ah3,'Position'), 'XAxisLocation','top',...
          'YTickLabel','','Color','none','XColor','k');
        bar(ah3,YZ_xpjt(:,1),YZ_xpjt(:,2),1,'g');
        xlabel('y distance (nm)'); xlim([roi_ymin,roi_ymax]);
      
        
      ah5 = subplot(2,5,[2,3]);
        pcolor(XZ_xx,XZ_zz,XZ); cbh=colorbar('East'); colormap jet; set(cbh,'Ycolor','r')
        rectangle('Position',[widc/2+px_sz,widc/2+px_sz,lenc,widc],'Curvature',1,'LineWidth',2,'LineStyle','--','EdgeColor','w');
        axis equal; axis([roi_xmin,roi_xmax,roi_zmin,roi_zmax]);
        set(ah5,'XTickLabel','','YTickLabel',''); %freezeColors;
      ah6 = subplot(2,5,[7,8]);
        nXZ_xpjt_x=(XZ_xpjt(:,1)-((widcpx/2+1)*px_sz))/(lenc/2);
        bar(ah6,nXZ_xpjt_x,XZ_xpjt(:,2),1,'g');
        nlim_min = (roi_xmin-(widc/2+px_sz))/(lenc/2);
        nlim_max = (roi_xmax-(widc/2+px_sz))/(lenc/2);
        xlim([nlim_min,nlim_max]); xlabel('norm position'); ylabel('Counts');
        hold on;
        ax6t =  axes('Position',get(ah6,'Position'), 'XAxisLocation','top',...
          'YTickLabel','','Color','none','XColor','k');
        bar(ah6,XZ_xpjt(:,1),XZ_xpjt(:,2),1,'g');
        xlabel('x distance (nm)'); xlim([roi_xmin,roi_xmax]);     
      
      ah8 = subplot(2,5,[4,5]);
        pcolor(XY_xx,XY_yy,XY); cbh=colorbar('East'); colormap jet; set(cbh,'Ycolor','r')
        rectangle('Position',[widc/2+px_sz,widc/2+px_sz,lenc,widc],'Curvature',1,'LineWidth',2,'LineStyle','--','EdgeColor','w');
        axis equal; axis([roi_xmin,roi_xmax,roi_ymin,roi_ymax]);
        set(ah8,'XTickLabel','','YTickLabel',''); %freezeColors;
      ah9 = subplot(2,5,[9,10]);
        nXY_xpjt_x=(XY_xpjt(:,1)-((widcpx/2+1)*px_sz))/(lenc/2);
        bar(ah9,nXY_xpjt_x,XY_xpjt(:,2),1,'g');
        nlim_min = (roi_xmin-(widc/2+px_sz))/(lenc/2);
        nlim_max = (roi_xmax-(widc/2+px_sz))/(lenc/2);
        xlim([nlim_min,nlim_max]); xlabel('norm position'); ylabel('Counts');
        hold on;
        ax9t =  axes('Position',get(ah9,'Position'), 'XAxisLocation','top',...
          'YTickLabel','','Color','none','XColor','k');
        bar(ah9,XY_xpjt(:,1),XY_xpjt(:,2),1,'g');
        xlabel('x distance (nm)'); xlim([roi_xmin,roi_xmax]);
      
        uicontrol('Style', 'text','String', currentdic,'Units','normalized','Position', [0.01 0.01 1 0.03]); 
      if savefig==1;
      saveppt2(pptfilname,'figure',fig_RI2,'scale',true,'stretch',false,'editmenu');
      end
      
%%  GENERATE SLICE ALONG X AXIS (YZ SLICE)
          fig_yz_RI = figure;
          ha = tight_subplot(2,ceil(bina/2),[.01 .03],[.1 .01],[.01 .01]);
          haaa=ceil(bina/2); habb=bina/2;
          if haaa ~= habb; delete(subplot(2,ceil(bina/2),bina+1));end;
      for yz_i=1:bina;
          % yz slice
          axes(ha(yz_i));  
          yzidx = find(roi_xminpx+(yz_i-1)*rimg_xstep <= rimg_x & rimg_x <= roi_xminpx+yz_i*rimg_xstep);
          rimg_yz_slice = [rimg_y(yzidx),rimg_z(yzidx)];
          RI_yz_slice = iQ_ReconImage(rimg_yz_slice,[roi_yranpx+1,roi_zranpx+1],mask_size,sigma,precision);
          
          imshow(RI_yz_slice,[]); set(gca,'ydir','normal'); colormap hot; title(['yz int hist,no. ',num2str(yz_i)]);
          rectangle('Position',[(widcpx/2+1)*10,(widcpx/2+1)*10,widcpx*10,widcpx*10],'Curvature',[1,1],'LineWidth',2,'LineStyle','--','EdgeColor','w');
          axis([roi_yminpx*10,roi_ymaxpx*10,roi_zminpx*10,roi_zmaxpx*10]);
          
          PF_yz{ii,yz_i}=rimg_yz_slice;
      end
      
          uicontrol('Style', 'text','String', currentdic,'Units','normalized','Position', [0.01 0.01 1 0.03]); 
          if savefig==1;
          saveppt2(pptfilname1,'figure',fig_yz_RI,'scale',true,'stretch',false,'editmenu','title',['Cell no. ' num2str(ii),'Protn Den:' num2str(ProtnDen) ' um-3']);
          end
          
      for yz_i=1:bina;  
          yzidx = find(roi_xminpx+(yz_i-1)*rimg_xstep <= rimg_x & rimg_x <= roi_xminpx+yz_i*rimg_xstep);
          xd = rimg_y(yzidx);
          yd = rimg_z(yzidx);
          xi = linspace(1,roi_yranpx+1,binb);
          yi = linspace(1,roi_zranpx+1,binc);
          xr = interp1(xi,0.5:numel(xi)-0.5,xd,'nearest');
          yr = interp1(yi,0.5:numel(yi)-0.5,yd,'nearest');
          
          Z = accumarray([yr xr]+0.5, 1, [binb binc]);    
          
          [yy, zz]=meshgrid(roi_yranpx/binb*px_sz:roi_yranpx/binb*px_sz:roi_yranpx*px_sz,roi_zranpx/binc*px_sz:roi_zranpx/binc*px_sz:roi_zranpx*px_sz);
          Z_xpjt = sum(Z); Z_ypjt = sum(Z,2);
          Z_xpjt_ax = yy(1,:)+0.5*yy(1,1); Z_ypjt_ax = zz(:,1)+0.5*zz(1,1); 
          
          fig_yzslice = figure;
          ah1 = subplot(4, 4, [1:3,5:7,9:11]);
          pcolor(yy,zz,Z); h=colorbar('East'); set(h,'Ycolor','r')
          title(['yz intensity histogram,no. ',num2str(yz_i)]); ylabel('z distance (nm)');
          rectangle('Position',[widc/2+px_sz,widc/2+px_sz,widc,widc],'Curvature',[1,1],'LineWidth',2,'LineStyle','--','EdgeColor','w');
          axis equal
          
          ah2 = subplot(4, 4, [13,14,15]);
          bar(ah2,Z_xpjt_ax,Z_xpjt,1,'g');
          xlabel('y distance (nm)'); ylabel('Counts');
          
          linkaxes([ah1 ah2 ],'x');
          pos2 = get(ah2,'Position');
          pos1 = get(ah1,'Position');
          pos2(3) = pos1(3);
          set(ah2,'Position',pos2)
          set(ah1,'XTickLabel','')
          pos2(2) = pos1(2) - pos2(4);
          set(ah2,'Position',pos2)
          
          ah3 = subplot(4, 4, [4,8,12]);
          barh(ah3,Z_ypjt_ax,Z_ypjt,1,'g'); 
          set (ah3, 'YAxisLocation', 'right', 'XAxisLocation', 'top','Color', 'w');
          
          linkaxes([ah1 ah3 ],'y');
          pos3 = get(ah3,'Position');
          pos3(4) = pos1(4);
          set(ah3,'Position',pos3)
          set(ah3,'YTickLabel','')
          pos3(1) = pos1(1) + pos1(3);
          set(ah3,'Position',pos3)
          
          ah4 = subplot(4,4,16);
          hold on; 
          rectangle('Position',[1+(yz_i-1)*rimg_xstep,1,rimg_xstep,2*widcpx],'LineWidth',1,'FaceColor','y')
          plot(rimg_x,rimg_y,'x','MarkerSize',2);
          axis([roi_xminpx,roi_xmaxpx,roi_yminpx,roi_ymaxpx]);
          set(ah4,'XTickLabel','','YTickLabel','','FontSize',6); xlabel('xy scatter plot');
          hold off; 
          rectangle('Position',[(widcpx/2+1),(widcpx/2+1),lencpx,widcpx],'Curvature',1,'LineWidth',2,'LineStyle','--','EdgeColor','b');
          
          pos4 = get(ah4,'Position');
          pos4(1) = pos3(1); pos4(2) = pos2(2); pos4(3) = pos3(3); pos4(4)=pos2(4);
          set(ah4,'Position',pos4)
          axis equal; 
          
          uicontrol('Style', 'text','String', currentdic,'Units','normalized','Position', [0.01 0.01 1 0.03]); 
          if savefig==1;
          saveppt2(pptfilname2,'figure',fig_yzslice,'scale',true,'stretch',false,'editmenu','title',['Cell no. ' num2str(ii),'Protn Den:' num2str(ProtnDen) ' um-3']);
          end
      end
      
%%  GENERATE SLICE ALONG Y AXIS (XZ SLICE)     
          fig_xz_RI = figure;
          
          ha = tight_subplot(ceil(bina/2),2,[.01 .03],[.1 .01],[.01 .01]);
          haaa=ceil(bina/2); habb=bina/2;
          if haaa ~= habb; delete(subplot(ceil(bina/2),2,bina+1));end;
      for xz_i=1:bina;
          % xz slice
          axes(ha(xz_i));        
          xzidx = find(roi_yminpx+(xz_i-1)*rimg_ystep <= rimg_y & rimg_y <= roi_yminpx+xz_i*rimg_ystep);
          rimg_xz_slice = [rimg_x(xzidx),rimg_z(xzidx)];
          RI_xz_slice = iQ_ReconImage(rimg_xz_slice,[roi_xranpx+1,roi_zranpx+1],mask_size,sigma,precision);
          
          imshow(RI_xz_slice,[]); set(gca,'ydir','normal'); colormap hot; title(['xz int hist,no. ',num2str(xz_i)]);
          rectangle('Position',[(widcpx/2+1)*10,(widcpx/2+1)*10,lencpx*10,widcpx*10],'Curvature',1,'LineWidth',2,'LineStyle','--','EdgeColor','w');
          axis([roi_xminpx*10,roi_xmaxpx*10,roi_zminpx*10,roi_zmaxpx*10]);
          
          PF_xz{ii,xz_i}=rimg_xz_slice;
      end
          uicontrol('Style', 'text','String', currentdic,'Units','normalized','Position', [0.01 0.01 1 0.03]); 
          if savefig==1;
          saveppt2(pptfilname1,'figure',fig_xz_RI,'scale',true,'stretch',false,'editmenu','title',['Cell no. ' num2str(ii),'Protn Den:' num2str(ProtnDen) ' um-3']);
          end
          
      for xz_i=1:bina;
          % xz slice
          xzidx = roi_yminpx+(xz_i-1)*rimg_ystep <= rimg_y & rimg_y <= roi_yminpx+xz_i*rimg_ystep;
          xd = rimg_x(xzidx);
          yd = rimg_z(xzidx);
          xi = linspace(1,roi_xranpx+1,binb);
          yi = linspace(1,roi_zranpx+1,binc);
          xr = interp1(xi,0.5:numel(xi)-0.5,xd,'nearest');
          yr = interp1(yi,0.5:numel(yi)-0.5,yd,'nearest');
                    
          Z = accumarray([yr xr]+0.5, 1, [binb binc]);    
          
          [xx, zz]=meshgrid(roi_xranpx/binb*px_sz:roi_xranpx/binb*px_sz:roi_xranpx*px_sz,roi_zranpx/binc*px_sz:roi_zranpx/binc*px_sz:roi_zranpx*px_sz);
          Z_xpjt = sum(Z); Z_ypjt = sum(Z,2);
          Z_xpjt_ax = xx(1,:)+0.5*xx(1,1); Z_ypjt_ax = zz(:,1)+0.5*zz(1,1); 
          
          fig_xzslice = figure;
          ah1 = subplot(4, 4, [1:3,5:7,9:11]);
          pcolor(xx,zz,Z); h=colorbar('East'); set(h,'Ycolor','r')
          title(['xz intensity histogram,no. ',num2str(xz_i)]); ylabel('z distance (nm)');
          rectangle('Position',[widc/2+px_sz,widc/2+px_sz,lenc,widc],'Curvature',1,'LineWidth',2,'LineStyle','--','EdgeColor','w');
          axis equal; 
          
          ah2 = subplot(4, 4, [13,14,15]);
          bar(ah2,Z_xpjt_ax,Z_xpjt,1,'g');
          xlabel('x distance (nm)'); ylabel('Counts');
          
          linkaxes([ah1 ah2 ],'x');
          pos2 = get(ah2,'Position');
          pos1 = get(ah1,'Position');
          pos2(3) = pos1(3);
          set(ah2,'Position',pos2)
          set(ah1,'XTickLabel','')
          pos2(2) = pos1(2) - pos2(4);
          set(ah2,'Position',pos2)
          
          ah3 = subplot(4, 4, [4,8,12]);
          barh(ah3,Z_ypjt_ax,Z_ypjt,1,'g'); 
          set (ah3, 'YAxisLocation', 'right', 'XAxisLocation', 'top','Color', 'w');
          
          linkaxes([ah1 ah3 ],'y');
          pos3 = get(ah3,'Position');
          pos3(4) = pos1(4);
          set(ah3,'Position',pos3)
          set(ah3,'YTickLabel','')
          pos3(1) = pos1(1) + pos1(3);
          set(ah3,'Position',pos3)
          
          ah4 = subplot(4,4,16);
          hold on; 
          rectangle('Position',[1,1+(xz_i-1)*rimg_ystep,lencpx+widcpx,rimg_ystep],'LineWidth',1,'FaceColor','y')
          plot(rimg_x,rimg_y,'x','MarkerSize',2);
          axis([roi_xminpx,roi_xmaxpx,roi_yminpx,roi_ymaxpx]); 
          set(ah4,'XTickLabel','','YTickLabel','','FontSize',6); xlabel('xy scatter plot');
          hold off;  
          rectangle('Position',[(widcpx/2+1),(widcpx/2+1),lencpx,widcpx],'Curvature',1,'LineWidth',2,'LineStyle','--','EdgeColor','b');
               
          pos4 = get(ah4,'Position');
          pos4(1) = pos3(1); pos4(2) = pos2(2); pos4(3) = pos3(3); pos4(4)=pos2(4);
          set(ah4,'Position',pos4)
          axis equal; 
          
          uicontrol('Style', 'text','String', currentdic,'Units','normalized','Position', [0.01 0.01 1 0.03]); 
          if savefig==1;
          saveppt2(pptfilname2,'figure',fig_xzslice,'scale',true,'stretch',false,'editmenu','title',['Cell no. ' num2str(ii),'Protn Den:' num2str(ProtnDen) ' um-3']);
          end
      end
 %%  GENERATE SLICE ALONG Z AXIS (XY SLICE)      
          fig_xy_RI=figure; 
          ha = tight_subplot(ceil(bina/2),2,[.01 .03],[.1 .01],[.01 .01]);
          haaa=ceil(bina/2); habb=bina/2;
          if haaa ~= habb; delete(subplot(ceil(bina/2),2,bina+1));end;
          
      for xy_i=1:bina;
          % xy slice
          axes(ha(xy_i));  
          xyidx = find(roi_zminpx+(xy_i-1)*rimg_zstep <= rimg_z & rimg_z <= roi_zminpx+xy_i*rimg_zstep);
          rimg_xy_slice = [rimg_x(xyidx),rimg_y(xyidx)];
          RI_xy_slice = iQ_ReconImage(rimg_xy_slice,[roi_xranpx+1,roi_yranpx+1],mask_size,sigma,precision);
          
          imshow(RI_xy_slice,[]); set(gca,'ydir','normal'); colormap hot; title(['xy int hist,no. ',num2str(xy_i)]);
          rectangle('Position',[(widcpx/2+1)*10,(widcpx/2+1)*10,lencpx*10,widcpx*10],'Curvature',1,'LineWidth',2,'LineStyle','--','EdgeColor','w');
          axis([roi_xminpx*10,roi_xmaxpx*10,roi_yminpx*10,roi_ymaxpx*10]);
          
          PF_xy{ii,xy_i}=rimg_xy_slice;
      end
          uicontrol('Style', 'text','String', currentdic,'Units','normalized','Position', [0.01 0.01 1 0.03]); 
          if savefig==1;
          saveppt2(pptfilname1,'figure',fig_xy_RI,'scale',true,'stretch',false,'editmenu','title',['Cell no. ' num2str(ii),'Protn Den:' num2str(ProtnDen) ' um-3']);
          end

      for xy_i=1:bina;
          % xy slice
          xyidx = roi_zminpx+(xy_i-1)*rimg_zstep <= rimg_z & rimg_z <= roi_zminpx+xy_i*rimg_zstep;
          xd = rimg_x(xyidx);
          yd = rimg_y(xyidx);
          xi = linspace(1,roi_xranpx+1,binb);
          yi = linspace(1,roi_yranpx+1,binc);
          xr = interp1(xi,0.5:numel(xi)-0.5,xd,'nearest');
          yr = interp1(yi,0.5:numel(xi)-0.5,yd,'nearest');
          
          Z = accumarray([yr xr]+0.5, 1, [binb binc]);    
          
          [xx, yy]=meshgrid(roi_xranpx/binb*px_sz:roi_xranpx/binb*px_sz:roi_xranpx*px_sz,roi_yranpx/binc*px_sz:roi_yranpx/binc*px_sz:roi_yranpx*px_sz);
          Z_xpjt = sum(Z); Z_ypjt = sum(Z,2);
          Z_xpjt_ax = xx(1,:)+0.5*xx(1,1); Z_ypjt_ax = yy(:,1)+0.5*yy(1,1); 
            
          fig_xyslice = figure;
          ah1 = subplot(4, 4, [1:3,5:7,9:11]);
          pcolor(xx,yy,Z); h=colorbar('East'); set(h,'Ycolor','r')
          title(['xy intensity histogram,no. ',num2str(xy_i)]); ylabel('y distance (nm)');
          rectangle('Position',[widc/2+px_sz,widc/2+px_sz,lenc,widc],'Curvature',1,'LineWidth',2,'LineStyle','--','EdgeColor','w');
          axis equal;
          
          ah2 = subplot(4, 4, [13,14,15]);
          bar(ah2,Z_xpjt_ax,Z_xpjt,1,'g');
          xlabel('x distance (nm)'); ylabel('Counts');
          
          linkaxes([ah1 ah2 ],'x');
          pos2 = get(ah2,'Position');
          pos1 = get(ah1,'Position');
          pos2(3) = pos1(3);
          set(ah2,'Position',pos2)
          set(ah1,'XTickLabel','')
          pos2(2) = pos1(2) - pos2(4);
          set(ah2,'Position',pos2)
          
          ah3 = subplot(4, 4, [4,8,12]);
          barh(ah3,Z_ypjt_ax,Z_ypjt,1,'g'); 
          set (ah3, 'YAxisLocation', 'right', 'XAxisLocation', 'top','Color', 'w');
          
          linkaxes([ah1 ah3 ],'y');
          pos3 = get(ah3,'Position');
          pos3(4) = pos1(4);
          set(ah3,'Position',pos3)
          set(ah3,'YTickLabel','')
          pos3(1) = pos1(1) + pos1(3);
          set(ah3,'Position',pos3)
          
          ah4 = subplot(4,4,16);
          hold on; 
          rectangle('Position',[1,1+(xy_i-1)*rimg_zstep,lencpx+widcpx,rimg_zstep],'LineWidth',1,'FaceColor','y')
          plot(rimg_x,rimg_z,'x','MarkerSize',2);
          axis([roi_xminpx,roi_xmaxpx,roi_zminpx,roi_zmaxpx]); 
          set(ah4,'XTickLabel','','YTickLabel','','FontSize',6); xlabel('xz scatter plot');  
          hold off;  
          rectangle('Position',[(widcpx/2+1),(widcpx/2+1),lencpx,widcpx],'Curvature',1,'LineWidth',2,'LineStyle','--','EdgeColor','b');
                    
          pos4 = get(ah4,'Position');
          pos4(1) = pos3(1); pos4(2) = pos2(2); pos4(3) = pos3(3); pos4(4)=pos2(4);
          set(ah4,'Position',pos4)
          axis equal; 
          
          uicontrol('Style', 'text','String', currentdic,'Units','normalized','Position', [0.01 0.01 1 0.03]); 
          if savefig==1;
          saveppt2(pptfilname2,'figure',fig_xyslice,'scale',true,'stretch',false,'editmenu','title',['Cell no. ' num2str(ii),'Protn Den:' num2str(ProtnDen) ' um-3']);
          end
      end
      end
  end
%% nest routine saveppt2
end