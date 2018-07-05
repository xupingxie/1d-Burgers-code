function gnodes=Gauss_1D_nodes_local(elem,n_gp)
% 
 a=elem(1);
 b=elem(2);
N=n_gp-1;
N1=N+1;
N2=N+2;

x=linspace(a,b,N1)';

y=cos((2*(0:N)'+1)*pi/(2*N+2))+(0.27/N1)*sin(pi*x*N/N2);
%y=cos(pi*((0:N)'-1/4)/(N1+0.5));
L=zeros(N1,N2);
Lp=zeros(N1,N2);

y0=1;
while(max(abs(y-y0))>eps)
    L(:,1)=1;
    Lp(:,1)=0;
    L(:,2)=y;
    Lp(:,2)=1;
    
    for k=2:N1
        L(:,k+1)=((2*k-1)*y.*L(:,k)-(k-1)*L(:,k-1))/k;
    end
    Lp=(N2)*(L(:,N1)-y.*L(:,N2))./(1-y.^2);
    y0=y;
    y=y0-L(:,N2)./Lp;
end
% % rule=n_gp;
% %   r = zeros(n_gp,1);
% % %  w = zeros(n_gp,1);
% % 
% %   if rule == 1
% %     r(1) = 0;
% % %    w(1) = 2;
% %   elseif rule == 2
% %     r(1) =-1.0d0 / sqrt(3.0d0);
% %     r(2) =-r(1);
% % %    w(1) = 1.0;
% % %    w(2) = 1.0;
% %   elseif rule == 3
% %     r(1) =-sqrt(3.0d0/5.0d0);
% %     r(2) = 0.0;
% %     r(3) =-r(1);
% % %    w(1) = 5.0d0 / 9.0d0;
% % %    w(2) = 8.0d0 / 9.0d0;
% % %    w(3) = w(1);
% %   elseif rule == 4
% %     r(1) =-sqrt((3.0d0+2.0*sqrt(6.0d0/5.0d0))/7.0d0);
% %     r(2) =-sqrt((3.0d0-2.0*sqrt(6.0d0/5.0d0))/7.0d0);
% %     r(3) =-r(2);
% %     r(4) =-r(1);
% % %    w(1) = 0.5d0 - 1.0d0 / ( 6.0d0 * sqrt(6.0d0/5.0d0) );
% % %    w(2) = 0.5d0 + 1.0d0 / ( 6.0d0 * sqrt(6.0d0/5.0d0) );
% % %    w(3) = w(2);
% % %    w(4) = w(1);
% %   elseif rule == 5
% %     r(1) =-sqrt(5.0d0+4.0d0*sqrt(5.0d0/14.0d0)) / 3.0d0;
% %     r(2) =-sqrt(5.0d0-4.0d0*sqrt(5.0d0/14.0d0)) / 3.0d0;
% %     r(3) = 0.0d0;
% %     r(4) =-r(2);
% %     r(5) =-r(1);
% % %    w(1) = 161.0d0/450.0d0-13.0d0/(180.d0*sqrt(5.0d0/14.0d0));
% % %    w(2) = 161.0d0/450.0d0+13.0d0/(180.d0*sqrt(5.0d0/14.0d0));
% % %    w(3) = 128.0d0/225.0d0;
% % %    w(4) = w(2);
% % %    w(5) = w(1);
% %   elseif rule == 6
% %     r(1) = -0.2386191861;
% %     r(2) = -0.6612093865;
% %     r(3) = -0.9324695142;
% %     r(4) = - r(1);
% %     r(5) = - r(2);
% %     r(6) = - r(3);
% % %    w(1) = .4679139346;
% % %     w(2) = .3607615730;
% % %     w(3) = .1713244924;
% % %     w(4) = w(1);
% % %     w(5) = w(2);
% % %     w(6) = w(3);
% %   else
% %     error('Quadrature rule not supported')
% %     keyboard
% %   end
% %   y=r;
  
 gnodes=(a*(1-y)+b*(1+y))/2;
  
end
