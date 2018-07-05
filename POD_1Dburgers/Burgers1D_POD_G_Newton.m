function [a0]=Burgers1D_POD_G_Newton(pod_u,FEM,C0,Mr,Sr,Ten,Nonlin0,Nonlin1,...
         f0,f1,Ones,opt,n_gp);
     
     converged=0;
     diverged = 0;
     iter=0;
     max_iter=20;
     tol_step=1.d-6;
     tol_resid=1.d-10;
     
     r=size(C0,1);
     
     a0=C0;
     
     Mat_A=Mr;
     Mat_B=opt.nu*Sr;
     
     A_f=Mat_A+opt.theta*opt.del*Mat_B;
     rhs_f=(-Mat_A+(1-opt.theta)*opt.del*Mat_B)*C0-...
         opt.del*(opt.theta*f1+(1-opt.theta)*f0);
     
     while(~converged && ~diverged)
         
         iter=iter+1;
         d_nonlin=sparse(r,r);
         Nh=sparse(r,1);
         
         if opt.theta~=0 || iter==1
             for i=1:r
                 Nonlin0(i,(i-1)*r+1:i*r)=(Ten(:,:,i)*a0)';
                 Nonlin1(i,(i-1)*r+1:i*r)=a0'*Ten(:,:,i);
             end
         end
         
         d_nonlin=d_nonlin+(Nonlin0+Nonlin1)*Ones;
         Nh=Nh+Nonlin0*Ones*a0;
         
         if iter==1
             rhs_f=rhs_f+(1-opt.theta)*opt.del*Nh;
         end
         
         Ac=A_f+opt.del*opt.theta*d_nonlin;
         bc=-(Mat_A*a0+opt.del*opt.theta*Mat_B*a0+opt.del*opt.theta*Nh+rhs_f);
         
         Newton_step=Ac\bc;
         C1=a0+Newton_step;
         
         norm_step=norm(Newton_step,2);
         norm_x=norm(C1,2);
         norm_resid=norm(Ac*Newton_step-bc,2);
         if iter==1
             norm_resid0=max(norm_resid,1e-20);
         end
         a0=C1;
         
         converged=(norm_step < tol_step) || (norm_resid< tol_resid && ...
             norm_step/norm_x<tol_step);
         diverged = (iter>max_iter);
         
         if diverged == 1
             fprintf(1,'Newton solver diverges\n');
         end
     end

end