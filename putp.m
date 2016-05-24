function ansputp=putp(arg1)
putpi=0;
ansputp=[0 0];
ax=axis;
xspacing=(ax(1,2)-ax(1,1))/100;
yspacing=(ax(1,4)-ax(1,3))/100;
hold on
while 1,
   putpi=putpi+1;
   [x,y]=ginput(1);
   
   if (x<ax(1)); break;
   end
   plot(x,y,'rx');
   number=num2str(putpi);
   text(x+xspacing,y+yspacing,number,'color','w')
   if nargin==1,   
   if arg1==0,
    Text=[' ',num2str(x),', ',num2str(y)];
   elseif arg1==1,
    Text=[' ',num2str(x)];
   elseif arg1==2, 
    Text=[' ',num2str(y)];
   end
   text(x,y,Text);
   end
   ansputp(putpi,:)=[x y];
end
hold off
end