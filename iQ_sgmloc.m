function [ XYWODC,XYWDC,XZWODC,XZWDC,YZWODC,YZWDC ] = iQ_sgmloc( DFT_pos,DFT_dr_ave )
% This function calculates the localization precision
% INPUTS:
% DFT_pos: position of markers
% DFT_dr_ave: average drift of position markers
% OUTPUTS: 
% XYWODC/XYWDC: localization histogram in XY without/with drift correction
% XZWODC/XZWDC: localization histogram in XZ without/with drift correction
% YZWODC/YZWDC: localization histogram in YZ without/with drift correction


%% PLOT OUT THE IMAGE W/ AND W/O DC (JUST FOR VISUAL PURPOSE)

[~, col]=size(DFT_pos); % find out how many reference spots  
XYWODC=cell(col/16,1); XZWODC=cell(col/16,1); YZWODC=cell(col/16,1);
XYWDC=cell(col/16,2);  XZWDC=cell(col/16,2);  YZWDC=cell(col/16,2);
go3D=abs(max(DFT_pos(:,5)));
 for i=1:col/16
     
     % remove outliers by intensity filter to take care of the Au blinking case 
     DFT_pos_raw_int = DFT_pos(:,16*(i-1)+15);
     intmean=mean(DFT_pos_raw_int);
     intstd=std(DFT_pos_raw_int);
     intul=intmean+4*intstd; intll=intmean-4*intstd;
     intidx=find(intll<DFT_pos_raw_int & DFT_pos_raw_int<intul);
     
     % x,y positions without drift correction
     xwodc=DFT_pos(intidx,16*(i-1)+3);%-min(DFT_pos(intidx,16*(i-1)+3));      
     ywodc=DFT_pos(intidx,16*(i-1)+4);%-min(DFT_pos(intidx,16*(i-1)+4));
     
     % x,y positions with drift correction
     xwdc=xwodc-DFT_dr_ave(intidx,3);      
     ywdc=ywodc-DFT_dr_ave(intidx,4);
     
     [x_xywodc,y_xywodc,M_xywodc]=iQ_locprec(xwodc,ywodc,20);    
     
     [x_xywdc,y_xywdc,M_xywdc]=iQ_locprec(xwdc,ywdc,20);
     results_xywdc=iQ_autoGaussianSurfML(x_xywdc,y_xywdc,M_xywdc);
      XYWODC{i,1}=[x_xywodc,y_xywodc,M_xywodc]; 
      XYWDC{i,1}=[x_xywdc,y_xywdc,M_xywdc]; XYWDC{i,2}=results_xywdc;
     figure;
     subplot(231)
     surf(x_xywodc,y_xywodc,M_xywodc); view(2);
     axis tight; 
     xlabel('x (nm)'); ylabel('y (nm)'); 
     title('xy before drift corr')
      subplot(234)
     surf(x_xywdc,y_xywdc,M_xywdc); view(2);
     axis tight;
     xlabel('x (nm)'); ylabel('y (nm)'); 
     title({'xy after drift corr';...
         ['\sigmax=',num2str(round(results_xywdc.sigmax*100)/100),'nm',...
         '\sigmay=',num2str(round(results_xywdc.sigmay*100)/100),'nm']});
     
     % if you have 3D data
     if go3D > 0
     
     zwodc=DFT_pos(intidx,16*(i-1)+5);%-min(DFT_pos(intidx,16*(i-1)+5));
     zwdc=zwodc-DFT_dr_ave(intidx,5);
     [x_xzwodc,z_xzwodc,M_xzwodc]=iQ_locprec(xwodc,zwodc,20);
     [y_yzwodc,z_yzwodc,M_yzwodc]=iQ_locprec(ywodc,zwodc,20);
     
     [x_xzwdc,z_xzwdc,M_xzwdc]=iQ_locprec(xwdc,zwdc,20);
     results_xzwdc=iQ_autoGaussianSurfML(x_xzwdc,z_xzwdc,M_xzwdc);
     [y_yzwdc,z_yzwdc,M_yzwdc]=iQ_locprec(ywdc,zwdc,20);
     results_yzwdc=iQ_autoGaussianSurfML(y_yzwdc,z_yzwdc,M_yzwdc);  
     XZWODC{i,1}=[x_xzwodc,z_xzwodc,M_xzwodc]; 
     YZWODC{i,1}=[y_yzwodc,z_yzwodc,M_yzwodc];     
     XZWDC{i,1}=[x_xzwdc,z_xzwdc,M_xzwdc]; XZWDC{i,2}=results_xzwdc;
     YZWDC{i,1}=[y_yzwdc,z_yzwdc,M_yzwdc]; YZWDC{i,2}=results_yzwdc;
          
     subplot(232)
     surf(x_xzwodc,z_xzwodc,M_xzwodc); view(2);
     axis tight; 
     xlabel('x (nm)'); ylabel('z (nm)'); 
     title('xz before drift corr')
     subplot(233)
     surf(y_yzwodc,z_yzwodc,M_yzwodc); view(2);
     axis tight; 
     xlabel('y (nm)'); ylabel('z (nm)'); 
     title('yz before drift corr')
     subplot(235)
     surf(x_xzwdc,z_xzwdc,M_xzwdc); view(2);
     axis tight;
     xlabel('x (nm)'); ylabel('z (nm)'); 
     title({'xz after drift corr';...
         ['\sigmax=',num2str(round(results_xzwdc.sigmax*100)/100),'nm',...
         '\sigmaz=',num2str(round(results_xzwdc.sigmay*100)/100),'nm']});
     subplot(236)
     surf(y_yzwdc,z_yzwdc,M_yzwdc); view(2);
     axis tight;
     xlabel('y (nm)'); ylabel('z (nm)'); 
     title({'yz after drift corr';...
         ['\sigmay=',num2str(round(results_yzwdc.sigmax*100)/100),'nm',...
         '\sigmaz=',num2str(round(results_yzwdc.sigmay*100)/100),'nm']});
     end
end
%% nested function iQ_locprec
function [xx,yy,M]=iQ_locprec(xxx,yyy,n)
%code from make 2d histogram
%n=20;  %number of bins in each dimension
subx=xxx; suby=yyy;
len=length(subx);

M=zeros(n,n);  %empty matrix that will become the 2d histogram

minx = min(subx); maxx = max(subx); rangex=maxx-minx;
stepx=rangex/n;

miny= min(suby); maxy=max(suby); rangey=maxy-miny;
stepy=rangey/n;

for j=1:len
    for x=1:n
        for y=1:n
            if subx(j)<=minx+x*stepx && subx(j)> minx+(x-1)*stepx && suby(j)<=miny+y*stepy && suby(j)>miny+(y-1)*stepy
                M(n-y+1,x)=M(n-y+1,x)+1;
            end
        end
    end
end

[xx, yy]=meshgrid(((minx+0.5*stepx):stepx:(maxx-0.5*stepx)),((maxy-0.5*stepy):-stepy:(miny+0.5*stepy)));
end
end

