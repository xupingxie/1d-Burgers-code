function [Ten,Ten_x0,Ten_0x]=POD_tensor_assemble_1D(r,pod_u,FEM,n_gp)

for i=1:r
    [T_x0, T_0x, T]=POD_tensor_trilinear_1D(pod_u,FEM,n_gp,i);    
	Ten(:,:,i)   = pod_u'*T*pod_u;
	Ten_x0(:,:,i) = pod_u'*T_x0*pod_u; 
	Ten_0x(:,:,i) = pod_u'*T_0x*pod_u;    
end

return

end
