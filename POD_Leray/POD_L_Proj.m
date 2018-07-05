%----------L2,Ritz projection
%----
format long
time_start=tic;
x_l=0;
x_r=1;
t=[0,1];

N=10000;
n=1024;


opt.nu=0.001;
opt.theta=0;
opt.del=1e-4;      
Tr=20;          % # of pod basis,  r
dr=10;       % # of projection basis, r1
podopt.L2_space=1;
podopt.plot_basis=0;

podopt.L2_proj=1;  % 1---L2 projection, 0.5--H1 semi ritz projection,0 H1 full proj

degree=1;
n_gp=5;



load u_fem_data.mat
 
 
 u0=u_fem(:,1);
 
% mesh=mesh_generator_1D(x_l,x_r,n);
% GDOF=global_dof_1D_fe_Lagrange(mesh,degree);
% FEM=struct('mesh',mesh,'GDOF',GDOF,'degree',degree);
% 
% 
% [Mass,Stiff]=FE_matrix_1D_Lagrange(a_name,FEM,n_gp);       %-----FE Mass, Stiff matrix

EH1=[];
EL2=[];
EHs=[];
Eer=[];
E=[];

p=0;
%for j=1:length(Tr)

podopt.p=Tr;

[pod_u,CumEng,CumEng_ratio,sum_lambda,lambda_p]=POD_basis(FEM,u_fem,Mass,Stiff,podopt);             %-------generate POD basis;

r=podopt.p;

Ten=POD_tensor_assemble_1D(r,pod_u,FEM,n_gp);        %---get the tensor;

Nonlin0 = zeros(r,r^2);
Nonlin1 = zeros(r,r^2);
Ones = zeros(r^2,1);
C=sparse(r,r);

for i=1:r
    Ones((i-1)*r+1:i*r,1:r)=eye(r);
end

Mr=pod_u'*Mass*pod_u;
Sr=pod_u'*Stiff*pod_u;

for j=1:length(dr);

C0=Mr\pod_u'*Mass*u0;

C(:,1)=C0;

save_index=100;


C_temp=zeros(dr(j),1);


M_pr=Mr(1:dr(j),:);
M_bar=Mr(1:dr(j),1:dr(j));
S_pr=Sr(1:dr(j),:);
S_bar=Sr(1:dr(j),1:dr(j));

if podopt.L2_proj==1
    Fv=M_pr;
    Mc=M_bar;
elseif podopt.L2_proj==0.5
    Fv=S_pr;
    Mc=S_bar;
else
    Fv=M_pr+S_pr;
    Mc=M_bar+S_bar;
end


for k=1:N
    
    C_bar=zeros(r,1);
    
    C_temp=Mc\(Fv*C0);
    
    C_bar(1:dr(j),1)= C_temp;
    
    
    Nh=zeros(r,1);
    
%    for i=1:r                      %--------my format
%        Nonlin0(i,(i-1)*r+1:i*r)=(Ten(:,:,i)*C_bar)';
%    end   
%    Nh=Nh+Nonlin0*Ones*C0;
 
    for i=1:r                   %----------zHU'S format
        Nh(i)=(C_bar)'*Ten(:,:,i)*C0;          %(u_bar u_x,\phi_i)
    end
    
    F=Mr*C0-opt.nu*opt.del*(Sr*C0)-opt.del*Nh;       %---F vector
    
    C1=Mr\F;
    if rem(k,save_index)==0
        C(:,k/save_index+1)=C1;
    end
    C0=C1;
end

p=p+1,
u_rom=pod_u*C;

Av=POD_AverDNS_new(u_rom,Mass);
E=[E,Av];   
%-----Caculate the error btw  FEM and POD ROM
% %     
 Err_L2=zeros(101,1);
% % %Err_H1=Err_L2;
% % %Err_H1_semi=Err_L2;
% % %Err_max=Err_L2;
% % 
for i=1:101
    
    [Err_max(i),Err_L2(i),Err_H1_semi(i),Err_H1(i)]=...
        POD_error_1D(Mass,Stiff,u_fem(:,i),u_rom(:,i));
end
% % 
 Err_L2_avg=sqrt(sum(Err_L2.^2)/101);
% % %Err_H1_avg=sqrt(sum(Err_H1.^2)/101);
% % %Err_H1_semi_avg=sqrt(sum(Err_H1_semi.^2)/101);
% % 
% % 
EL2=[EL2,Err_L2_avg];
%    EH1=[EH1,Err_H1_avg];
    %Eer=[Eer,CumEng_ratio(r)];
%end    
end

% Zmax= max(max(u_fem))*1.5;
% Zmin= min(min(u_fem))*1.5;
% t_initial=0;
% t_final=1;
% dt=100/N;
% t_checkpts=0:dt:1;
% x=GDOF.P_g;
% figure(j)
% mesh(x,t_checkpts', u_rom');
% xlabel('x'); ylabel('t');
% view([-1 -1 1]);
% axis([t_initial,t_final,x(1),x(end),Zmin,Zmax]);
% title(sprintf('L-ROM-projection r=%d,R=%d',Tr,dr),'Fontsize',20);

ddd=POD_AverDNS_new(u_fem,Mass);
cputime=toc(time_start)

Zmax= max(max(u_fem))*1.5;
Zmin= min(min(u_fem))*1.5;
t_initial=0;
t_final=1;

dt=100/N;
t_checkpts=0:dt:1;
x=GDOF.P_g;
figure(1)
mesh(x,t_checkpts', u_rom');
xlabel('x'); ylabel('t');
view([-1 -1 1]);
axis([t_initial,t_final,x(1),x(end),Zmin,Zmax]);
title(['L-Proj r=', num2str(Tr), ' and r1= ', num2str(dr)],'Fontsize',15);
%colorbar
% r=zeros(length(E));
% r=Av;
% figure(2)
% subplot(2,1,1)
% plot(dr,E,'bo--','Linewidth',2.5)
% hold on
% plot([0 max(dr)],[ddd ddd],'r-','Linewidth',2.5)
% hold off
% legend('Leray-DF','DNS');
% title(['Average energy vs delta at r= ', num2str(Tr)])
% xlabel(['$\delta','$'],'Interpreter','latex','Fontsize',15);
% ylabel('Average energy','Fontsize',15);
% subplot(2,1,2)
% plot(dr,EL2,'bo--','Linewidth',2.5)
% title(['Average errors vs delta at r= ', num2str(Tr)])

