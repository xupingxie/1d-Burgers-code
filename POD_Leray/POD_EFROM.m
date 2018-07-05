%----------Differential fitler process
%----
format long
clear all;
time_start=tic;
x_l=0;
x_r=1;
t=[0,1];

N=10000;
n=1024;

opt.nu=0.001;

opt.del=1e-4;    %----Delta t 

%Tr=2:4:50;      % # of pod basis

Tr=[2,4,6,8,10,12,14];     %----number of POD basis,  r =10;

podopt.L2_space=1;  % 1 L2, 0 H1 basis
podopt.plot_basis=0;
EL2=[];
EH1=[];
EHs=[];
save_index=100;

degree=1;
n_gp=5;         %---number of quadrature points


load u_fem_data.mat


u0=u_fem(:,1);    %---initial value u0


podopt.p=Tr;

[pod_u,CumEng,CumEng_ratio,POD_all,Diag_S,lp]=POD_basis(FEM,u_fem,Mass,Stiff,podopt);             %-------generate POD basis;   

r=podopt.p;

[~,Ten,~]=POD_tensor_assemble_1D(r,pod_u,FEM,n_gp);        %---get the tensor;

Nonlin0 = zeros(r,r^2);
Nonlin1 = zeros(r,r^2);
Ones = zeros(r^2,1);
C=sparse(r,r);

for i=1:r
    Ones((i-1)*r+1:i*r,1:r)=eye(r);
end

Mr=pod_u'*Mass*pod_u;

Sr=pod_u'*Stiff*pod_u;

C0=Mr\pod_u'*Mass*u0;

C_temp=C0;
C(:,1)=C_temp;
%C0=C(:,1);

%de=0.000001;     %----filter radius, delta
%de=0:0.00004:0.002;  %-----optimal delta=0.000012;
%de=0.00052; %r=10 &15
de=0.00069;
%de=0.000375;
%de=0.0007;
%de=0.00181:0.00002:0.00189;
%de=0.0006;
%de=0.002;
%de=0;
p=0;
E=[];   %------store the aver data

for j=1:length(Tr)
    
    C0=C(:,1);
    
    Err_L2=sparse(101,1),

    for k=1:N
        
        Nh=zeros(r,1);
%         for i=1:r
%             Nonlin0(i,(i-1)*r+1:i*r)=(Ten(:,:,i)*C0)';
%         end
%         
%         Nh=Nh+Nonlin0*Ones*C0;
        for i=1:r
            Nh(i) = (C0)'*Ten(:,:,i)*C0; %
        end        

        F=Mr*C0-opt.nu*opt.del*(Sr*C0)-opt.del*Nh;       %---F vector
        
        C1_tmp=Mr\F;                        %--forward Euler, compute w^n+1
        C1=(Mr+de(j)^2*Sr)\(Mr*C1_tmp);         %----DF, compute w_bar ^n+1
        
        if rem(k,save_index)==0
            C(:,k/save_index+1)=C1;
        end
        C0=C1;
    end
    
    u_rom=pod_u*C;

%plot(GDOF.P_g,u_rom(:,end),'ob');


%-----Caculate the error btw  FEM and POD ROM
Av=POD_AverDNS_new(u_rom,Mass),
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
p=p+1,
    
end

ddd=POD_AverDNS_new(u_fem,Mass);
cputime=toc(time_start),

Zmax= max(max(u_rom))*1.5;
Zmin= min(min(u_rom))*1.5;
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
%title(['EF-ROM with r=20 $\delta=0.000375','$'],'Interpreter','latex','Fontsize',15);
title(['EF-ROM r=', num2str(Tr)],'Fontsize',15);
% %colorbar
% % r=zeros(length(E));
% % r=Av;
% figure(2)
% plot(de,E,'b','Linewidth',2.5)
% hold on
% plot([0 0.002],[20.3218 20.3218],'r-','Linewidth',2.5)
% hold off
% legend('EF-ROM-DF','DNS');
% title('Optimal delta EF-ROM-DF r=20','Fontsize',15)
% xlabel('R1','Fontsize',15);
% ylabel('Average KE','Fontsize',15);
% %xlabel(['$\delta','$'],'Interpreter','latex','Fontsize',15);
% % ylabel(['$\frac{1}{101}\sum_{j=1}^{101}\|u_j\|_{L^2}','$'],'Interpreter','latex','Fontsize',15);
% % %ylabel('L2-Average','Fontsize',15);
% % 

% figure(1)
% plot(de,E,'b','Linewidth',2.5)
% hold on
% plot([0 max(de)],[ddd ddd],'r','Linewidth',2.5)
% hold off
% legend('EF-ROM','DNS');
% title('Average KE vs delta with r=20','Fontsize',15);
% %xlabel('r1','Fontsize',15);
% xlabel(['$\delta','$'],'Interpreter','latex','Fontsize',15);
% ylabel(['$\frac{1}{101}\sum_{j=1}^{101}\|u_j\|_{L^2}','$'],'Interpreter','latex','Fontsize',15);
% 
% figure(2)
% plot(de,EL2,'bo--','Linewidth',2.5)
% title(['Average errors vs delta at r=20'],'Fontsize',15);
% %xlabel('r1','Fontsize',15);
% xlabel(['$\delta','$'],'Interpreter','latex','Fontsize',15);
% ylabel(['$\frac{1}{101}\sum_{j=1}^{101}\|e_j\|_{L^2}','$'],'Interpreter','latex','Fontsize',15);
% 
