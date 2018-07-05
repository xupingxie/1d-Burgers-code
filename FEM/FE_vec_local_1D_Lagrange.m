function v_l=FE_vec_local_1D_Lagrange(f_name,t,elem,degree,d_index,n_gp)
v_l=zeros(degree+1,1);

nodes=Gauss_1D_nodes_local(elem,n_gp);
weights=Gauss_1D_weights_local(elem,n_gp);
m=f_name(nodes);
%m=feval(f_name,nodes,t);   %---time-dependent;
for i=1:degree+1
    fl=m.*shape_fun_1D_Lagrange(nodes,elem,...
            degree,i,d_index);
    v_l(i)=sum(fl.*weights);
end
return

