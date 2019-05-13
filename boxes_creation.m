function [x_box,y_box,z_box,box_length,box_count_x,box_count_y,box_count_z]=boxes_creation(Lx,Ly,Lz,rho)

x=1.01:0.001:2;
spacing=x*rho;
box_count_x=Lx./spacing;
a=floor(box_count_x);
b=box_count_x./a;
[~,n]=min(b);
% n=1;

x=x(n);
spacing=x*rho;
box_count_x=round(Lx/spacing);
box_count_y=round(Ly/spacing);
box_count_z=round(Lz/spacing);

box_length=Lx/box_count_x;

x_box=linspace(0,Lx,box_count_x);
y_box=linspace(0,Lz,box_count_y);
z_box=linspace(0,Ly,box_count_z);