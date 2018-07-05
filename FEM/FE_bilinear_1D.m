function M_l=FE_bilinear_1D(ker,phi,test,w_g)

M_l= test'*diag(w_g.*ker)*phi;

return