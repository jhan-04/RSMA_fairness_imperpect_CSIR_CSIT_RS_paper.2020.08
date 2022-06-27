function [MA_p,tou_p, P1_p,P2_p, Pc_p,Rs_p]=RS_noma(gam_dB,rho,P)
MA_p=2;
P2_p=0;


    v=0;
    
    %%%%%%%%%%%%%channel
    gam=10^(gam_dB/20);
        
        

        theta=acos(1-2*rho);
        
        
        h1=1/sqrt(2)*[1;1];
        h2=(gam)/sqrt(2)*[1;exp(-1i*theta)];

        %%%%%%%%%%%%%%%%%%%%%%%
         v=0;
    
    
    for t=0:1/200:1
        v=v+1;
        h1_=h1/sqrt(1+norm(h1)^2*rho*P*t);
       % h1_=h1/sqrt(1+norm(h1)^2*P*t);
        h2_=h2;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        A=h1_*h1_';
        B=h2_*h2_';
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cvx_begin quiet
        variable F(2,2) hermitian
        variable x(1)
        minimize(-x)
        subject to
        trace(A*F)>=x;
        trace(B*F)>=x;
        trace(F)==1;
        F== semidefinite(2);
        cvx_end
        %%%%%%%%%%%%%%%%%%%%%%%%%%
trace(A*F);
trace(B*F);
        
        [U,W,Z] = svds(F);%V=U*W*Z'
        obj_max=-100;
        n=length(W);
        for m=1:1000
            r=sqrt(1/2)*(randn(n,1)+1i*randn(n,1));%sqrt(var/2)*(randn(1,N)+1i*randn(1,N))
            v_=U*sqrt(W)*r;
            v_=v_/sqrt(v_'*v_);
            obj=min(v_'*A*v_,v_'*B*v_);
            if obj>=obj_max
                obj_max=obj;
                v_max=v_;
            end
        end
        f_c=v_max;
        norm(f_c)^2;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        Rc=min(log2(1+abs(h1_'*f_c)^2*(1- t)*P),log2(1+abs(h2_'*f_c)^2*(1- t)*P));
        log2(1+abs(h1_'*f_c)*(1- t)*P);
        log2(1+abs(h2_'*f_c)*(1- t)*P);
    R(v)=  Rc+ log2(1+norm(h1)^2*rho*P*t);
    end
        
    t=0:1/200:1;

    Rs_1=max(R);
    k=find(Rs_1==R);
    t1=t(k);
   
    
    
    
    
    
    
    
    tou_p=t1;
    P2_p=0;
    P1_p=P*t1;
    Pc_p=(1-t1)*P;
    Rs_p=Rs_1;
end
    