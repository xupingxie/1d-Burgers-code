function KEplot(u_rom,u_fem,Tr,dr)

[m,n]=size(u_rom);

E=[];
E1=[];
for i = 1:n
    y=norm(u_rom(:,i),2);
    y1=norm(u_fem(:,i),2);
    E=[E,y];
    E1=[E1,y1];
end



p=1:1:101;

figure(1)
plot(p,E,'b--','Linewidth',2.5);
hold on
plot(p,E1,'r-*','Linewidth',2.5);
title(sprintf('EF-ROM-proj r=%d,R=%d',Tr,dr),'Fontsize',15);
legend('EF-ROM-proj','DNS');
%xlabel(['$\delta=0.5','$'],'Interpreter','latex','Fontsize',15);
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
