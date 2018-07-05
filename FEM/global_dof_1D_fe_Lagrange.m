function GDOF=global_dof_1D_fe_Lagrange(mesh,degree)

n=length(mesh.P);
N=degree;
P_g=[];

P_g=linspace(mesh.P(1),mesh.P(n),n+(n-1)*(N-1));

T_g(:,1)=1:N+1;

for i=2:n-1
    T_g(:,i)=T_g(N+1,i-1):T_g(N+1,i-1)+N;
end

GDOF=struct('P_g',P_g,'T_g',T_g);
return
end