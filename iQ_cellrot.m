function [ celltheta ] = iQ_cellrot( DIC, mask_Red_grid,dil_sz,IM_info,PC_sgmErflt )
%UNTITLED Summary of this function goes here
%INPUTS:
%   DIC: DIC image
%   mask_Red_grid: mask x,y coordinates generated from iQ_imgshift
%   dil_size: dilate size for mask dilation in iQ_cellfnd
%   IM_info: pixel information
%   PC_sgmErflt: PC data after filtering
%OUTPUTS:
%   celltheta:  a five column cell, 
%               column 1: rotation angle in degree
%               column 2: corner x, corner y,length and width of cell in nm
%               column 3: volume of the cell in um^3
%% TOP of the routine
px_x=IM_info(1); px_y=IM_info(2);

% ROTATE THE CELL IN XY PLANE, find out the rotation angles for every cells using cellmask
  Redmask=zeros(size(DIC));
  Redmask(sub2ind(size(Redmask),mask_Red_grid(:,2),mask_Red_grid(:,1)))=1;
  cc = bwconncomp(Redmask, 4);
  celltheta=cell(cc.NumObjects,3);
  for ii=1:cc.NumObjects
      close all
      PC_sgmErflt_mat=PC_sgmErflt{ii};
      if isempty(PC_sgmErflt_mat); % skip cells contains no fitting results
          continue
      else
      cellidx=cc.PixelIdxList{ii};
      cellorientation=false(size(Redmask));
      cellorientation(cellidx)=true;
      se = strel('disk',dil_sz);  
      cellorientation=imerode(cellorientation,se);
      [cell_y,cell_x]=find(cellorientation~=0);
    
      STATS = regionprops(cellorientation,'Centroid','MajorAxisLength','MinorAxisLength', 'Orientation');
      
      if length(STATS)~=1;
          STATS = STATS(1);
      end
      if STATS.Orientation >0; STATS.Orientation = 180- STATS.Orientation;
      end
      % assume the cell has the cylinder with 2 hemisphere cap shape
      cell_vol=(4/3+(STATS.MajorAxisLength*px_x-STATS.MinorAxisLength*px_y))*pi*(STATS.MinorAxisLength/2*px_y)^2*10^-9;
      
      figure; subplot(121);plot(cell_x,cell_y,'o'); 
      axis equal; grid;axis([min(cell_x)-15,max(cell_x)+15,min(cell_y)-15,max(cell_y)+15]); % add or minus 15 just for visual prupose
      title({sprintf('Centroid = (%0.1f,%0.1f); Orientation (in degree) = %0.1f', STATS.Centroid(1),STATS.Centroid(2),STATS.Orientation);sprintf('MaxAxis = %0.1f um ; MinAxis = %0.1f um ',STATS.MajorAxisLength*px_x/1000,STATS.MinorAxisLength*px_y/1000)});
      gca;
      line([0,gca],[STATS.Centroid(2),STATS.Centroid(2)]);
      
      m=1;
      while 1
      theta=input('rotae angle (in degree)= ');
      theta=pi/180*theta;
      %rotate the coordinates and convert back to cartesian coordinates
      
      % deal all rotation in nm instead of pixel
      r=((cell_x*px_x).^2+(cell_y*px_y).^2).^0.5;
      t=atan((cell_y*px_y)./(cell_x*px_x));
      t2=t-theta;
      cell_x_new=r.*cos(t2);
      cell_y_new=r.*sin(t2);
      %centroid_new = [sum(cell_x_new)/size(cell_x_new,1) sum(cell_y_new)/size(cell_y_new,1)];
      [rectx,recty,~,~] = minboundrect(cell_x_new,cell_y_new);
      cenx=mean(rectx(1:4)); ceny=mean(recty(1:4));     centroid_new=[cenx,ceny];
      lenfit=[sqrt((rectx(1)-rectx(2))^2+(recty(1)-recty(2))^2),sqrt((rectx(2)-rectx(3))^2+(recty(2)-recty(3))^2)];
      lenx=max(lenfit); leny=min(lenfit);
      
      subplot(122); plot(cell_x_new,cell_y_new,'o'); axis image;
      %rectangle('Position',[centroid_new(1)-STATS.MajorAxisLength/2*px_x,centroid_new(2)-STATS.MinorAxisLength/2*px_y, STATS.MajorAxisLength*px_x,STATS.MinorAxisLength*px_y],'Curvature',1,'LineWidth',2,'LineStyle','--')
       rectangle('Position',[centroid_new(1)-lenx/2,centroid_new(2)-leny/2, lenx,leny],'Curvature',1,'LineWidth',2,'LineStyle','--')
      title({sprintf('Centroid new = (%0.2f,%0.2f); Cell Volume = %0.2f um^3; ',centroid_new(1),centroid_new(2),cell_vol); sprintf('MaxAxis = %0.1f um ; MinAxis = %0.1f um ',lenx/1000,leny/1000)});
      line([0,gca],[centroid_new(2),centroid_new(2)]);
      axis equal; grid;axis([min(cell_x_new)-10*px_x,max(cell_x_new)+10*px_x,min(cell_y_new)-10*px_y,max(cell_y_new)+10*px_y]);
      
      a=input('is the cell flat (Y/N)? ','s');
      A=upper(a);
      if A=='N'; m=m+1; end;
      if A=='Y',break; end; 
      end
 
       celltheta{ii,1}=theta;  
       celltheta{ii,2}=[centroid_new(1)-lenx/2,centroid_new(2)-leny/2, lenx,leny];
       celltheta{ii,3}=cell_vol;
      end
  end
%% nested function minboundrect
end

