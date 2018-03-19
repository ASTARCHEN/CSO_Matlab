clc
clear
close all

fun = inline('sqrt(x(:,1).^2+x(:,2).^2)','x');
[res, swarm] = cso(fun,100,[-8,-8],[8,8],1000,0.2);
xi = -8:0.05:8;
yi = -8:0.05:8;
[x,y] = meshgrid(xi, yi);
z = sqrt(x.^2+y.^2);
mesh(x,y,z);
hold on
plot3(swarm(1),swarm(2), res, 'r+')