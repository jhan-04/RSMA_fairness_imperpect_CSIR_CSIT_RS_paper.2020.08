function [MA_x,t_x, P1_x,P2_x, Pc_x,Rs_x,fc_x]=multi_sdr_imperfect(P,h1,h2,rho,beta)


MA_x=4;
P2_x=0;
P1_x=0;
Pc_x=P;
t_x=0;
Nt=length(h1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A=h1*h1';
B=h2*h2';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if rho==0
    
    fc_x=h1/norm(h1);
    
else
    %%%%%%%%%%%%%%%%%%%SDR
    cvx_begin quiet
    variable F(Nt,Nt) hermitian
    variable x(1)
    minimize(-x)
    subject to
    trace(A*F)>=x;
    trace(B*F)>=x;
    trace(F)==1;
    F== hermitian_semidefinite(Nt);
    cvx_end
    %%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%randaomization
    
    [U,W,Z] = svds(F);%V=U*W*Z'
    obj_max=-100;
    n=length(W);
    for m=1:10000
        r=sqrt(1/2)*(randn(Nt,1)+1i*randn(Nt,1));%sqrt(var/2)*(randn(1,N)+1i*randn(1,N))
        v_=U*sqrt(W)*r;
        v_=v_/sqrt(v_'*v_);
        obj=min(v_'*A*v_,v_'*B*v_);
        if obj>=obj_max
            obj_max=obj;
            v_max=v_;
        end
    end
    fc_x=v_max;
    
    
    Rc1=log2(1+abs(h1'*fc_x)^2*Pc_x/(1+beta*P+norm(h1)^2*rho*P1_x));
    Rc2=log2(1+abs(h2'*fc_x)^2*Pc_x/(1+beta*P+norm(h2)^2*rho*P2_x));
    Rs_x=min(Rc1,Rc2);
    
end

end