
t=[0:1/5:360];
t=t.*pi/180;
% без поворота
xel0 = 0 + 257.*cos(t);
yel0 = 0 + 150.*sin(t);
% с поворотом
phi=0;
for i=1:4
     phi=phi+45; % угол поворота градусы
     xel = 0 + 257.*cos(t)*cosd(phi)+150.*sin(t)*sind(phi);
     yel = 0 + 150.*sin(t)*cosd(phi)-257.*cos(t)*sind(phi);
hold on
plot(xel,yel,'linewidth',2),grid on
 
end
plot(257.*cos(t),257.*sin(t),'linewidth',2)
axis equal