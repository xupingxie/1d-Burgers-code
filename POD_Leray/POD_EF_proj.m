%----------Differential fitler process
%----
format long
time_start = tic;
x_l=0;
x_r=1;
t=[0,1];

N=10000;
n=1024;

opt.nu=0.001;

opt.del=1e-4;    %----Delta t 

%Tr=2:4:50;      % # of pod basis

Tr=10;     %----number of POD basis,  r;
dr=10;   %-----projection basis, r1;

podopt.L2_space=1;  % 1 L2, 0 H1 basis
podopt.plot_basis=0;
EL2=[];
EH1=[];
EHs=[];
save_index=100;

degree=1;
n_gp=5;         %---number of quadrature points


load u_fem_data.mat

p=0;
E=[];   %------store the aver data

u0=u_fem(:,1);    %---initial value u0


podopt.p=Tr;

[pod_u,CumEng,CumEng_ratio,POD_all,Diag_S,lp]=POD_basis(FEM,u_fem,Mass,Stiff,podopt);             %-------generate POD basis;   

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

%for j=1:length(de)

C_temp=zeros(dr(j),1);


M_pr=Mr(1:dr(j),:);
M_bar=Mr(1:dr(j),1:dr(j));
S_pr=Sr(1:dr(j),:);
S_bar=Sr(1:dr(j),1:dr(j));

%if podopt.L2_proj==1
    Fv=M_pr;
    Mc=M_bar;


    C0=C(:,1);

    for k=1:N
        
%        Nh=sparse(r,1);
        Nh=zeros(r,1);
%        for i=1:r                         %-----my format
%            Nonlin0(i,(i-1)*r+1:i*r)=(Ten(:,:,i)*C0)';
%        end         
%        Nh=Nh+Nonlin0*Ones*C0;


         for i=1:r                   %----------zHU'S format
             Nh(i)=(C0)'*Ten(:,:,i)*C0;
         end

        F=Mr*C0-opt.nu*opt.del*(Sr*C0)-opt.del*Nh;       %---F vector
        
        C1_tmp=Mr\F;                        %--forward Euler, compute w^n+1
        
        
        C_temp=Mc\(Fv*C1_tmp);             %projection;
        C_bar=zeros(r,1);  
        C_bar(1:dr(j),1)= C_temp;
        
        
        if rem(k,save_index)==0
            C(:,k/save_index+1)=C_bar;
        end
        C0=C_bar;
    end
    
u_rom=pod_u*C;

%plot(GDOF.P_g,u_rom(:,end),'ob');
%KEplot(u_rom,u_fem,Tr,dr);
p=p+1,

Av=POD_AverDNS_new(u_rom,Mass);
E=[E,Av]; 

Err_L2=zeros(101,1);

% Err_H1=Err_L2;
% % Err_H1_semi=Err_L2;
% % Err_max=Err_L2;
% % 
for i=1:101
    
    [Err_max(i),Err_L2(i),Err_H1_semi(i),Err_H1(i)]=...
        POD_error_1D(Mass,Stiff,u_fem(:,i),u_rom(:,i));
end
% % 
Err_L2_avg=sqrt(sum(Err_L2.^2)/101),
% % Err_H1_avg=sqrt(sum(Err_H1.^2)/101);
% % Err_H1_semi_avg=sqrt(sum(Err_H1_semi.^2)/101);
% % 
EL2=[EL2,Err_L2_avg],
% % EH1=[EH1,Err_H1_avg];
% % EHs=[EHs,Err_H1_semi_avg];

end

ddd=POD_AverDNS_new(u_fem,Mass);
cputime=toc(time_start),
Zmax= max(max(u_fem))*1.5;
Zmin= min(min(u_fem))*1.5;
t_initial=0;
t_final=1;
dt=100/N;
t_checkpts=0:dt:1;
x=GDOF.P_g;
figure(j)
mesh(x,t_checkpts', u_rom');
xlabel('x'); ylabel('t');
view([-1 -1 1]);
axis([t_initial,t_final,x(1),x(end),Zmin,Zmax]);
title(sprintf('EF-ROM-projection r=%d,r_1=%d',Tr,dr(j)),'Fontsize',15);



% figure(2)
% plot(dr,E,'b','Linewidth',2.5)
% hold on
% plot([0 5],[20.3218 20.3218],'r-','Linewidth',2.5)
% hold off
% legend('L-Proj','DNS');
% title('Optimal r1 L-Proj r=6','Fontsize',15)
% xlabel('R1','Fontsize',15);
% ylabel('Average KE','Fontsize',15);
% figure(1)
% plot(dr,E,'b','Linewidth',2.5)
% hold on
% plot([0 max(dr)],[ddd ddd],'r','Linewidth',2.5)
% hold off
% legend('EF-Proj','DNS');
% title('Average KE vs r_1 with r=20','Fontsize',15);
% xlabel('r1','Fontsize',15);
% %xlabel(['$\delta','$'],'Interpreter','latex','Fontsize',15);
% ylabel(['$\frac{1}{101}\sum_{j=1}^{101}\|u_j\|_{L^2}','$'],'Interpreter','latex','Fontsize',15);
% figure(2)
% plot(dr,EL2,'bo--','Linewidth',2.5)
% title(['Average errors vs r_1 at r=20'],'Fontsize',15);
% xlabel('r1','Fontsize',15);
% %xlabel(['$\delta','$'],'Interpreter','latex','Fontsize',15);
% ylabel(['$\frac{1}{101}\sum_{j=1}^{101}\|e_j\|_{L^2}','$'],'Interpreter','latex','Fontsize',15);
% 

