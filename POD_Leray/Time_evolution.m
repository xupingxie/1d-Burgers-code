
%function E=Time_evolution(u_rom)

[m,n]=size(u_rom);

E_rom=[];
E_fem=[];
for i = 1:n
    %y=norm(u_rom(:,i),2);
    y=sqrt((u_rom(:,i))'*Mass*u_rom(:,i));
    %y1=norm(u_fem(:,i),2);
    y1 = sqrt((u_fem(:,i))'*Mass*u_fem(:,i));
    E_rom=[E_rom,y];
    E_fem=[E_fem,y1];
end

r=1:1:101;

figure(1)
plot(r,E_rom,'b--','Linewidth',2.5);
hold on
plot(r,E_fem,'r-*','Linewidth',2.5);
%title(['EF-ROM with r=20 $\delta=0.000375','$'],'Interpreter','latex','Fontsize',15);
title('Time evolution of EF-Proj with r=10, r_1=10','Fontsize',15);
legend('EF-Proj','DNS')
xlabel('time','Fontsize',15);
ylabel(['$\|u_j\|_{L^2}','$'],'Interpreter','latex','Fontsize',15);

% figure(2)
% plot(r,E1,'r--');


%plot(de,E,'b',de,r,'r--','Linewidth',2)
% plot(de,EL2,'b-*');
% title('Optimal delta EF-ROM','Fontsize',20);
% legend('EF-ROM','DNS')
% xlabel(['$\delta','$'],'Interpreter','latex','Fontsize',20);
% ylabel('L2-Error','Fontsize',20);
