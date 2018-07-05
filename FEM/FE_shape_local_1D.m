function [phi,phi_x,weights]=FE_shape_local_1D(elem,degree,n_gp)

nodes=Gauss_1D_nodes_local(elem,n_gp);
weights=Gauss_1D_weights_local(elem,n_gp);
N=length(nodes);
phi=zeros(N,degree+1);
phi_x=zeros(N,degree+1);
for i=1:degree+1
    phi(:,i)=shape_fun_1D_Lagrange(nodes,elem,degree,i,0);
    phi_x(:,i)=shape_fun_1D_Lagrange(nodes,elem,degree,i,1);
end
    
end