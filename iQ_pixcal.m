function [px_x,px_y]=iQ_pixcal(x_cal,y_cal,d_l)
%% THIS FUNCTION CALCULATES THE DISTANCE PERPIXEL IN X AND Y
% INPUTs:
% x_cal: DIC calibration image in x
% y_cal: DIC calibration image in y
% d_l: distance per line of calibration slide in um, ex: 40lp/mm means d_l=12.5um.
% OUTPUTs:
% px_x: nm/pixel for x
% px_y: nm/pixel for y

%% main
x_cal=double(x_cal); y_cal=double(y_cal);
[r_x,~]=size(x_cal); [r_y,~]=size(y_cal);

px_idx_x=(1:1:r_x)'; px_idx_y=(1:1:r_y)';
lx1=round(1*r_x/6);lx2=round(3*r_x/6);lx3=round(5*r_x/6);
ly1=round(1*r_y/6);ly2=round(3*r_y/6);ly3=round(5*r_y/6);
cx1=x_cal(lx1,:)';
cx2=x_cal(lx2,:)';
cx3=x_cal(lx3,:)';
cy1=y_cal(:,ly1);
cy2=y_cal(:,ly2);
cy3=y_cal(:,ly3);

x_ave=(cx1+cx2+cx3)/3; y_ave=(cy1+cy2+cy3)/3;
midx=px_idx_x; midx(:)=mean(x_ave);midy=px_idx_y; midy(:)=mean(y_ave);
P_x = InterX([px_idx_x';x_ave'],[px_idx_x';midx']);
P_y = InterX([px_idx_y';y_ave'],[px_idx_y';midy']);
dis_x=(0:d_l:(length(P_x(1,:))-1)*d_l);
dis_y=(0:d_l:(length(P_y(1,:))-1)*d_l);
polyfitx = polyfit(P_x(1,:),dis_x,1);
polyfity = polyfit(P_y(1,:),dis_y,1);
px_x=polyfitx(1)*1000; px_y=polyfity(1)*1000;
fitx_x=px_idx_x; fity_x=polyfitx(1)*fitx_x+polyfitx(2);
fitx_y=px_idx_y; fity_y=polyfity(1)*fitx_y+polyfity(2);

figure;
subplot(321); imshow(x_cal,[]);
line([1,r_x],[lx1,lx1]),line([1,r_x],[lx2,lx2]),line([1,r_x],[lx3,lx3])
title('x pixel calibration')
subplot(322);imshow(y_cal,[]);
line([ly1,ly1],[1,r_y]),line([ly2,ly2],[1,r_y]),line([ly3,ly3],[1,r_y])
title('y pixel calibration')
subplot(323); plot(px_idx_x,cx1,px_idx_x,cx2,px_idx_x,cx3);
legend('LineProfile 1','LineProfile 2','LineProfile 3');
xlabel('x (px)'); ylabel('Intenisty (a.u.)');
subplot(324); plot(px_idx_y,cy1,px_idx_y,cy2,px_idx_y,cy3);
legend('LineProfile 1','LineProfile 2','LineProfile 3');
xlabel('y (px)'); ylabel('Intenisty (a.u.)');
subplot(325); plot(P_x(1,:),dis_x,'o',fitx_x,fity_x,'r-'); grid;
legend('x pixel data','fit curve');
xlabel('x (px)'); ylabel('Distance (\mum)');
title(['x calibration: ',num2str(px_x), 'nm'])
subplot(326); plot(P_y(1,:),dis_y,'o',fitx_y,fity_y,'r-'); grid;
legend('y pixel data','fit curve');
xlabel('y (px)'); ylabel('Distance (\mum)');
title(['y calibration: ',num2str(px_y), 'nm'])

%% nested function
function P = InterX(L1,varargin)
%INTERX Intersection of curves
%   P = INTERX(L1,L2) returns the intersection points of two curves L1 
%   and L2. The curves L1,L2 can be either closed or open and are described
%   by two-row-matrices, where each row contains its x- and y- coordinates.
%   The intersection of groups of curves (e.g. contour lines, multiply 
%   connected regions etc) can also be computed by separating them with a
%   column of NaNs as for example
%
%         L  = [x11 x12 x13 ... NaN x21 x22 x23 ...;
%               y11 y12 y13 ... NaN y21 y22 y23 ...]
%
%   P has the same structure as L1 and L2, and its rows correspond to the
%   x- and y- coordinates of the intersection points of L1 and L2. If no
%   intersections are found, the returned P is empty.
%
%   P = INTERX(L1) returns the self-intersection points of L1. To keep
%   the code simple, the points at which the curve is tangent to itself are
%   not included. P = INTERX(L1,L1) returns all the points of the curve 
%   together with any self-intersection points.
%   
%   Example:
%       t = linspace(0,2*pi);
%       r1 = sin(4*t)+2;  x1 = r1.*cos(t); y1 = r1.*sin(t);
%       r2 = sin(8*t)+2;  x2 = r2.*cos(t); y2 = r2.*sin(t);
%       P = InterX([x1;y1],[x2;y2]);
%       plot(x1,y1,x2,y2,P(1,:),P(2,:),'ro')

%   Author : NS
%   Version: 3.0, 21 Sept. 2010

%   Two words about the algorithm: Most of the code is self-explanatory.
%   The only trick lies in the calculation of C1 and C2. To be brief, this
%   is essentially the two-dimensional analog of the condition that needs
%   to be satisfied by a function F(x) that has a zero in the interval
%   [a,b], namely
%           F(a)*F(b) <= 0
%   C1 and C2 exactly do this for each segment of curves 1 and 2
%   respectively. If this condition is satisfied simultaneously for two
%   segments then we know that they will cross at some point. 
%   Each factor of the 'C' arrays is essentially a matrix containing 
%   the numerators of the signed distances between points of one curve
%   and line segments of the other.

    %...Argument checks and assignment of L2
    error(nargchk(1,2,nargin));
    if nargin == 1,
        L2 = L1;    hF = @lt;   %...Avoid the inclusion of common points
    else
        L2 = varargin{1}; hF = @le;
    end
       
    %...Preliminary stuff
    x1  = L1(1,:)';  x2 = L2(1,:);
    y1  = L1(2,:)';  y2 = L2(2,:);
    dx1 = diff(x1); dy1 = diff(y1);
    dx2 = diff(x2); dy2 = diff(y2);
    
    %...Determine 'signed distances'   
    S1 = dx1.*y1(1:end-1) - dy1.*x1(1:end-1);
    S2 = dx2.*y2(1:end-1) - dy2.*x2(1:end-1);
    
    C1 = feval(hF,D(bsxfun(@times,dx1,y2)-bsxfun(@times,dy1,x2),S1),0);
    C2 = feval(hF,D((bsxfun(@times,y1,dx2)-bsxfun(@times,x1,dy2))',S2'),0)';

    %...Obtain the segments where an intersection is expected
    [i,j] = find(C1 & C2); 
    if isempty(i),P = zeros(2,0);return; end;
    
    %...Transpose and prepare for output
    i=i'; dx2=dx2'; dy2=dy2'; S2 = S2';
    L = dy2(j).*dx1(i) - dy1(i).*dx2(j);
    i = i(L~=0); j=j(L~=0); L=L(L~=0);  %...Avoid divisions by 0
    
    %...Solve system of eqs to get the common points
    P = unique([dx2(j).*S1(i) - dx1(i).*S2(j), ...
                dy2(j).*S1(i) - dy1(i).*S2(j)]./[L L],'rows')';
              
    function u = D(x,y)
        u = bsxfun(@minus,x(:,1:end-1),y).*bsxfun(@minus,x(:,2:end),y);
    end
end
end








