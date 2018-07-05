function mesh=mesh_generator_1D(xl,xr,n)
% mesh generator function
P=zeros(1,n+1);
T=zeros(2,n);
E=zeros(3,n+1);

for i=0:n
    P(i+1)=((n-i)*xl+(i*xr))/n;
end

for i=1:n
    T(1,i)=i;
    T(2,i)=i+1;
end
E(1,:)=1:n+1;
for i=1:n+1
    if (i==1)
        E(2,i)=0;
        E(3,i)=1;
    elseif (i==n+1)
        E(2,i)=n;
        E(3,i)=0;
    else
        E(2,i)=i-1;
        E(3,i)=i;
    end
end

mesh=struct('P',P,'T',T,'E',E);
return; 

