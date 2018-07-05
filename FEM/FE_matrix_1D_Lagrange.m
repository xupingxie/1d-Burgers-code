function [M,M_x0, M_0x, S]=FE_matrix_1D_Lagrange(a_name,FEM,n_gp)

mesh=FEM.mesh;
GDOF=FEM.GDOF;
degree=FEM.degree;

M=sparse(size(GDOF.P_g,2),size(GDOF.P_g,2));
M_0x=M;
M_x0=M;
S=sparse(size(GDOF.P_g,2),size(GDOF.P_g,2));

if a_name==1
    ker=ones(n_gp,1);
end

for k=1:size(mesh.T,2)
    
    elem=mesh.P(mesh.T(:,k));
    
    [phi,phi_x,wg]=FE_shape_local_1D(elem,degree,n_gp);
    
    M_l=FE_bilinear_1D(ker,phi,phi,wg);
	M_x0l=FE_bilinear_1D(ker,phi,phi_x,wg);
	M_0xl=FE_bilinear_1D(ker,phi_x,phi,wg);
    S_l=FE_bilinear_1D(ker,phi_x,phi_x,wg);

    M(GDOF.T_g(:,k),GDOF.T_g(:,k))=M(GDOF.T_g(:,k),GDOF.T_g(:,k))+M_l;
	M_0x(GDOF.T_g(:,k),GDOF.T_g(:,k))=M_0x(GDOF.T_g(:,k),GDOF.T_g(:,k))+M_0xl;
	M_x0(GDOF.T_g(:,k),GDOF.T_g(:,k))=M_x0(GDOF.T_g(:,k),GDOF.T_g(:,k))+M_x0l;    
	S(GDOF.T_g(:,k),GDOF.T_g(:,k))=S(GDOF.T_g(:,k),GDOF.T_g(:,k))+S_l;
end
return;

end
