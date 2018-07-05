%load u_fem_data.mat

function [Av_eng]=POD_AverDNS_new(u_fem,mass)

[m, n]=size(u_fem);

E1=[];

for i=1:n
%     e=norm(u_fem(:,i),2);
    %e = sqrt((u_fem(:,i))'*mass*u_fem(:,i));
	e = (u_fem(:,i))'*mass*u_fem(:,i);
    E1=[E1,e];    
end


%Sum_eng=sum(E1);
Av_eng=sum(E1)/n;

