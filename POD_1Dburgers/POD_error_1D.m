function [Err_max,Err_L2,Err_H1_semi,Err_H1]=POD_error_1D(M,S,u_rom,u_fem)

u_h=u_fem;
u_pod=u_rom;
% mesh=FEM.mesh;
% GDOF=FEM.GDOF;
% degree=FEM.degree;

% function [ Err_max, Err_L2, Err_H1_semi, Err_H1]= FEM_POD_error_MS...
%     (M, S, u_h, u_pod)
% calculate errors between finite element solution and the POD-ROM
% solution on the same mesh

Err_max        = max(abs(u_h-u_pod));
Err_L2_square      = (u_h-u_pod)'*M*(u_h-u_pod);
Err_H1_semi_square = (u_h-u_pod)'*S*(u_h-u_pod);
Err_H1_square      = Err_L2_square+Err_H1_semi_square;

Err_L2      = sqrt(Err_L2_square);
Err_H1_semi = sqrt(Err_H1_semi_square);
Err_H1      = sqrt(Err_H1_square);


end