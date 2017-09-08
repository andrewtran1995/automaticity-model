x=rand(10,5);
y=rand(10,5);
z=rand(10,5);
t=0:4;
h=uicontrol(gcf,'style','slider','units','pix','position',[100 5 300 20]);
set(h,'min',min(t),'max',max(t));
set(h,'callback','i=find(t==nearest(get(h,''value'')));plot3(x(:,i),y(:,i),z(:,i))');
for i=1:length(t)
plot3(x(:,i),y(:,i),z(:,i));
pause(1);
set(h,'value',t(i));
end