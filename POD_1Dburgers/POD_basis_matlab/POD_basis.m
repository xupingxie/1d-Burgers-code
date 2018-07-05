function [POD_u,CumEng,CumEng_ratio,POD_all,Diag_S,Lambda_sequence]=POD_basis(FEM,u_snap,Mass,Stiffness,podopt)

GDOF=FEM.GDOF;
%--- generate POD basis
% compute basis vectors by efficient SVD method

if podopt.L2_space==1
    M = Mass;     % basis in L2 space
else
    M = Mass + Stiffness; % basis in H1 space
end
[POD_all,Diag_S,d,CumEng,CumEng_ratio] = PODbasis_generator...
    (u_snap,M);
%save(PODfile, 'POD_all','Diag_S','d','CumEng','CumEng_ratio')

%%---calculate the phi_1norm;

% Mr=POD_all'*Mass*POD_all;
% Sr=POD_all'*S*POD_all;
% l2=diag(Mr);
% ls=diag(Sr);
% v1=l2+ls;
% ph1=sum(v1.*(Diag_S.^2));


%--- calcualte percentage
ratio_p    = CumEng_ratio(podopt.p);
% ratio_p: ( sum_{j=1}^r \lambda_j )/( sum_{j=1}^d \lambda_j )
sum_lambda = sum(Diag_S.^2);
% sum_lambda: sum_{j=1}^d \lambda_j
% Lambda_p: sum_{j=r+1}^d \lambda_j
Lambda_sequence = sqrt(sum_lambda-CumEng);

% fprintf(1,['Rank of snapshots d= ',num2str(d),'\n']);
% if d<podopt.p
%     fprintf(1,'Rank of snapshots d< p\n');
%     podopt.p = d;
% end
% 
% switch podopt.type
%     case 'fix'
%         r   = podopt.p;
%         Lambda_p   = sqrt(sum_lambda-C1(r));
%     case 'tol'
%         [~, r] = min(abs(Lambda_sequence-podopt.tol));
%         Lambda_p   = sqrt(sum_lambda-C1(r));
%     case 'seq'
%         %--- find a sequence of r such that lambda dacay by half starting from 'Lambda_p'
%         r = podopt.p;
%         Lambda_p   = sqrt(sum_lambda-C1(r));
%         r_sequence = find_r_sequence(podopt.numseq,Lambda_sequence,Lambda_p);
%     case 'Rtol' % ratio of kinetic energy
%         [~, r] = min(abs(1-C2-podopt.tol));
%         Lambda_p   = sqrt(sum_lambda-C1(r));
% end
r=podopt.p;
lambda_p=sqrt(sum_lambda-CumEng(r));

POD_u = POD_all(:,1:r);

% fprintf(1,['# of POD basis in ROM is ', num2str(r), '\n']);

if podopt.plot_basis==1
figure(1)   
    for i=1:4,
        aa= subplot(2,2,i);
        pl=2*i-1;
        plot(GDOF.P_g,POD_all(:,pl));
        title(['$\varphi','_',num2str(pl),'$'],'Interpreter','latex');
        h = get(aa, 'title');
        set(h, 'Fontsize', 16, 'FontWeight', 'bold')
    end
end

bound_H1 = 0;
bound_L2 = 0;
%bound_L2= error in L2, ||u-\sum_j=1^r(u,phi_j)*phi_j||
for k = r+1:d
    bound_H1 = bound_H1 + (POD_all(:,k)'*S*POD_all(:,k))*Diag_S(k).^2;
    bound_L2 = bound_L2 + (POD_all(:,k)'*Mass*POD_all(:,k))*Diag_S(k).^2;
end
bound_H1 = sqrt(bound_H1);
bound_L2 = sqrt(bound_L2);

end
%--- plot basis functions
% [nm,~] = size(Mass);
% [n ,m] = size(u_save);
% dim    = n/nm;
% 
% if podopt.plot_basis == 1 % plot first few POD modes
%     for k = 1: dim
%         POD_title = [' DQ= ',num2str(podopt.dq)];
%         fig_num = plot_POD_basis(x,e_conn,POD_x((k-1)*nm+1:k*nm,:),fig_num,POD_title);
%         Figure_PODmodes = [File_fig, num2str(k), '_basis'];
%         print('-f', '-depsc', Figure_PODmodes)
%     end
% end

