function y=shape_fun_1D_Lagrange(x,elem,degree,...
    shape_index,d_index)
%%------updated on Dec 3, 2014;

N=length(x);
k=shape_index;
a=elem(1);
b=elem(2);


r1=(x-b)/(a-b);
r2=(x-a)/(b-a);
dr1=1/(a-b);
dr2=1/(b-a);
if degree==1
    if d_index==0
        if k==1
            f=r1;   %%lambda1
        else
            f=r2;   %%lambda2
        end
    else
        if k==1
            f=dr1;
        else
            f=dr2;
        end
    end
elseif degree==2
    if d_index==0
        if k==1
            f=2*r1.^2-r1;
        elseif k==2
            f=4*r1.*r2;
        else
            f=2*r2.^2-r2;
        end
    else
        if k==1
            f=4*r1.*dr1-dr1;
        elseif k==2
            f=4*r2.*dr1+4*r1.*dr2;
        else
            f=4*r2.*dr2-dr2;
        end
    end
elseif degree==3
    if d_index==0
        if k==1
            f=(1/2)*r1.*(3*r1-2).*(3*r1-1);
        elseif k==2
            f=(9/2)*r1.*r2.*(3*r1-1);
        elseif k==3
            f=(9/2)*r1.*r2.*(3*r2-1);
        else
            f=(1/2)*r2.*(3*r2-2).*(3*r2-1);
        end
    elseif d_index==1
        if k==1
            f=(1/2)*(27*r1.^2.*dr1-18*r1.*dr1+2*dr1);
        elseif k==2
            f=(9/2)*((dr1.*r2++r1.*dr2).*(3*r1-1)+3*r1.*r2.*dr1);
        elseif k==3
            f=(9/2)*((dr1.*r2++r1.*dr2).*(3*r2-1)+3*r1.*r2.*dr2);
        else
            f=(1/2)*(27*r2.^2.*dr2-18*r2.*dr2+2*dr2);
        end
    end
end

y=f;

end
