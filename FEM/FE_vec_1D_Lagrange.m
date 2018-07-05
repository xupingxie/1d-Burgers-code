function v=FE_vec_1D_Lagrange(f_name,t,FEM,d_index,n_gp)

%mesh=mesh_generator_1D(xl,xr,n);
mesh=FEM.mesh;
GDOF=FEM.GDOF;
degree=FEM.degree;
%GDOF=global_dof_1D_fe_Lagrange(mesh,degree);

v=zeros(size(GDOF.P_g,2),1);
for k=1:size(GDOF.T_g,2)
    elem=mesh.P(mesh.T(:,k));
    v_l=FE_vec_local_1D_Lagrange(f_name,t,elem,degree,...
        d_index,n_gp);
    v(GDOF.T_g(:,k))=v(GDOF.T_g(:,k))+v_l;
end
return;

end
