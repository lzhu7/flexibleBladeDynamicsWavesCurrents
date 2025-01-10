function [J]=Jacob(F_Sys,X)

%input arguments:
%F_Sys--system of nonlinear functions
%X--computing points
%output:
%J--Jacobian matrix

Nx=12;
Dx=1e-6;
J=zeros(6,12);%6 equations 12 variables
%forward finite differences
for j=1:Nx
    Xj=zeros(Nx,1);
    Xj(j)=Dx;
    J(:,j)=(F_Sys(X+Xj)-F_Sys(X))/Dx;
end