function [T_x0,T_0x,T]=POD_tensor_trilinear_1D(pod_u,FEM,n_gp,flag)

mesh=FEM.mesh;
GDOF=FEM.GDOF;
degree=FEM.degree;

p=flag;

T=sparse(size(GDOF.P_g,2),size(GDOF.P_g,2));
T_x0=T;
T_0x=T;

for k=1:size(mesh.T,2)
    
    elem=mesh.P(mesh.T(:,k));
    
    u_local(:,1)=GDOF.T_g(:,k);
    
    [phi,phi_x,wg]=FE_shape_local_1D(elem,degree,n_gp);
    
    ker=phi*pod_u(u_local,flag);
    Tx0_loc = FE_bilinear_1D(ker,phi_x,phi,wg);   %----(ker*h_x,h)= (phi(i) h_x,h) -> (phi(i) phi_x,phi) -> (u_x u, phi(i)) = (uu_x, phi(i));

	T_loc  = FE_bilinear_1D(ker,phi,phi,wg);

	%ker1=phi_x*pod_u(u_local,flag);               %-----(ker1*h,h)-> (phi_x(i) h,h)-> (phi_x(i) phi, phi)-> (uu,phi_x(i)) 
	T0x_loc = FE_bilinear_1D(ker,phi,phi_x,wg);
    
    T(GDOF.T_g(:,k),GDOF.T_g(:,k))   = T(GDOF.T_g(:,k),GDOF.T_g(:,k))+T_loc;
	T_x0(GDOF.T_g(:,k),GDOF.T_g(:,k)) = T_x0(GDOF.T_g(:,k),GDOF.T_g(:,k))+Tx0_loc;
	T_0x(GDOF.T_g(:,k),GDOF.T_g(:,k)) = T_0x(GDOF.T_g(:,k),GDOF.T_g(:,k))+T0x_loc;
    
end

return

end

    
