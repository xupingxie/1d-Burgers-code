
format long
x_l=0;
x_r=1;
t=[0,1];

N=100;
save_index=1;
n=1024;


opt.nu=0.001;
opt.theta=1;
opt.del=1e-2;
podopt.p=10;      % # of pod basis
podopt.L2_space=1;
r=podopt.p;
podopt.plot_basis=0;
record_evstep=0;
record_evsnap=1-record_evstep;

degree=1;
n_gp=5;


f_name='f_fun';
u0_name='u0_fun';
true_u_name='true_u_fun';
true_dx_name='true_dx_fun';
a_name=1;
u_dbc_name='u_dbc_fun';


load u_fem_BkE.mat

u_snap=u_fem;
u0=u_snap(:,1);

mesh=mesh_generator_1D(x_l,x_r,n);
GDOF=global_dof_1D_fe_Lagrange(mesh,degree);
FEM=struct('mesh',mesh,'GDOF',GDOF,'degree',degree);


[Mass,Stiff]=FE_matrix_1D_Lagrange(a_name,FEM,n_gp);       %-----FE Mass, Stiff matrix

[pod_u,CumEng]=POD_basis(FEM,u_snap,Mass,Stiff,podopt);             %-------generate POD basis;

r=podopt.p;
Ten=POD_tensor_assemble_1D(r,pod_u,FEM,n_gp);        %---get the tensor;

Nonlin0 = zeros(r,r^2);
Nonlin1 = zeros(r,r^2);
Ones = zeros(r^2,1);


for i=1:r
    Ones((i-1)*r+1:i*r,1:r)=eye(r);
end

Mr=pod_u'*Mass*pod_u;
Sr=pod_u'*Stiff*pod_u;

C0=Mr\pod_u'*Mass*u0;

nt=1;
C(:,nt)=C0;


% t=0;d_index=0;
% f0=pod_u'*FE_vec_1D_Lagrange(f_name,t,FEM,d_index,n_gp);
f0=0;
p=0;


 for k=1:N
     
     %t=opt.del*k;
     %f1=pod_u'*FE_vec_1D_Lagrange(f_name,t,FEM,0,n_gp);
     f1=0;
     
     [c1]=Burgers1D_POD_G_Newton(pod_u,FEM,C0,Mr,Sr,Ten,Nonlin0,Nonlin1,...
         f0,f1,Ones,opt,n_gp);
%      
%      if record_evstep==1
%          nt=nt+1;
%          C(:,nt)=c1;
%      elseif record_evsnap==1
%          nt=nt+1;
%          C(:,nt)=c1;
%      end
     C0=c1;
     if rem(k,save_index)==0
         %v_fe=u_new;        
         %eval(['save v_',int2str(k),'.txt v_fe -ascii'])
         C(:,k/save_index+1)=c1;
     end
    p=p+1,
 end

u_rom=pod_u*C;

plot(GDOF.P_g,u_rom(:,end),'ob')



Err_L2=zeros(101,1);
Err_H1=Err_L2;
Err_H1_semi=Err_L2;
Err_max=Err_L2;

for i=1:101
    
    [Err_max(i),Err_L2(i),Err_H1_semi(i),Err_H1(i)]=...
        POD_error_1D(Mass,Stiff,u_fem(:,i),u_rom(:,i));
end

Err_L2_avg=sqrt(sum(Err_L2.^2))/101,
Err_H1_avg=sqrt(sum(Err_H1.^2))/101,


