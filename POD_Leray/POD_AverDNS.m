%load u_fem_data.mat

function Av=POD_AverDNS(u_fem,C)

[m, n]=size(C);

E=[];

for i=1:n
    e=sum(C(:,i).^2);
    E=[E,e];    
end

Av=sum(E)/n;
%sqrt(u_fem'*C*u_fem);
%Av=sum(E)/n;

end

