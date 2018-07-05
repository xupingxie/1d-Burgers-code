%----------Differential fitler process
%----
format long
x_l=0;
x_r=1;
t=[0,1];

N=10000;
n=1024;

opt.nu=0.001;
opt.theta=0;
opt.del=1e-4;
%Tr=2:4:50;      % # of pod basis
Tr=6;
podopt.L2_space=1;  % 1 L2, 0 H1 basis
podopt.plot_basis=0;
EL2=[];
EH1=[];
EHs=[];
save_index=100;

degree=1;
n_gp=5;

load u_fem_data.mat

u0=u_fem(:,1);

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

C0=Mr\pod_u'*Mass*u0;

C_temp=C0;

 
C(:,1)=C_temp;
C0=C(:,1);

%de=0:0.002:0.2; %optimal de = 0.135
%de=0.02:0.02:1;
%de=0.2:0.004:0.4;
%de=0.4:0.04:0.8;
%de=1:0.004:1.4;
%de=0.00001:0.00002:0.002;  %-----optimal delta=0.000012;
%de=0.00181:0.00002:0.00189;
%de=0.0002
%de=0.14;
de = 0.01;
p=0;
E=[];
for j=1:length(de)
    
    C0=C(:,1);
    for k=1:N
        
        C_bar=(Mr+(de(j)^2)*Sr)\(Mr*C0);
       
        Nh=sparse(r,1);
        
%         for i=1:r                     %----old format
%             Nonlin0(i,(i-1)*r+1:i*r)=(Ten(:,:,i)*C_bar)';
%         end        
%         Nh=Nh+Nonlin0*Ones*C0;

 
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
    
   u_rom=pod_u*C;
%   plot(GDOF.P_g,u_rom(:,end),'ob');
    
   Av=POD_AverDNS(u_rom,C),
   E=[E,Av];
    
    %-----Caculate the error btw  FEM and POD ROM
    
     Err_L2=zeros(101,1);
% %     Err_H1=Err_L2;
% %     Err_H1_semi=Err_L2;
% %     Err_max=Err_L2;
% %     
    for i=1:101
        
        [Err_max(i),Err_L2(i),Err_H1_semi(i),Err_H1(i)]=...
            POD_error_1D(Mass,Stiff,u_fem(:,i),u_rom(:,i));
    end
% % %     
     Err_L2_avg=sqrt(sum(Err_L2.^2)/101);
% % %     Err_H1_avg=sqrt(sum(Err_H1.^2)/101);
% % %     Err_H1_semi_avg=sqrt(sum(Err_H1_semi.^2)/101);
% % %     
      EL2=[EL2,Err_L2_avg];
%     EH1=[EH1,Err_H1_avg];
%     EHs=[EHs,Err_H1_semi_avg];
    p=p+1,
    
end

Zmax= max(max(u_fem))*1.5;
Zmin= min(min(u_fem))*1.5;
t_initial=0;
t_final=1;

% dt=100/N;
% t_checkpts=0:dt:1;
% x=GDOF.P_g;
% figure(1)
% mesh(x,t_checkpts', u_rom');
% xlabel('x'); ylabel('t');
% view([-1 -1 1]);
% axis([t_initial,t_final,x(1),x(end),Zmin,Zmax]);
% title(['L-DF with r=10 $\delta=0.16','$'],'Interpreter','latex','Fontsize',15);
%plot(GDOF.P_g,POD_u(:,))
%colorbar
% r=zeros(length(E));
% r=Av;
% figure(1)
% plot(de,E,'b','Linewidth',2.5)
% hold on
% plot([0 0.2],[20.3218 20.3218],'r-','Linewidth',2.5)
% hold off
% legend('Leray-DF','DNS');
% title('Optimal delta L-DF with r=6','Fontsize',15)
% xlabel(['$\delta','$'],'Interpreter','latex','Fontsize',15);
% ylabel(['$\|u_j\|_{L^2}','$'],'Interpreter','latex','Fontsize',15);
% % 

% figure(1)
% plot(de,E,'b--','Linewidth',2.5);
% hold on
% plot(r,E1,'r-*','Linewidth',2.5);
% title(['EFROM-DF with r=10 $\delta=0.002','$'],'Interpreter','latex','Fontsize',15);
% title('KE evolution of L-proj with r=20,r1=18','Fontsize',15);
% legend('L-proj','DNS')
% xlabel('time','Fontsize',15);
% ylabel(['$\|u_j\|_{L^2}','$'],'Interpreter','latex','Fontsize',15);
