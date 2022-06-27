
function  [fc]=Find_fc(P,p,h1,h2,rho,N0, sigma_e1,sigma_e2)




Nt=length(h1);

h1_=h1/sqrt(norm(h1)^2*rho*p(1)+sigma_e1^2*P+N0);
h2_=h2/sqrt(norm(h2)^2*rho*p(2)+sigma_e2^2*P+N0);
A=h1_*h1_';
B=h2_*h2_';


if rho==0
    
    fc=h1/norm(h1);
    
else
    %%%%%%%%%%%%%%%%%%%SDR
    cvx_begin quiet
    variable F(Nt,Nt) hermitian
    variable x(1)
    minimize(-x)
    subject to
    trace(A*F)>=x;
    trace(B*F)>=x;
    trace(F)<=1;
    F== semidefinite(Nt);
    cvx_end
    %%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%randaomization
    
    [U,W,Z] = svds(F);%V=U*W*Z'
    obj_max=-100;
    n=length(W);
    for m=1:1000000
        r=sqrt(1/2)*(randn(Nt,1)+1i*randn(Nt,1));%sqrt(var/2)*(randn(1,N)+1i*randn(1,N))
        v_=U*sqrt(W)*r;
        v_=v_/sqrt(v_'*v_);
        obj=min(v_'*A*v_,v_'*B*v_);
        if obj>=obj_max
            obj_max=obj;
            v_max=v_;
        end
    end
    fc=v_max;
    
% 
% if rho==0
%     
%     fc=h1/norm(h1);
% else
% h1_=h1/sqrt(norm(h1)^2*rho*p(1)+sigma_e1^2*P+N0);
% h2_=h2/sqrt(norm(h2)^2*rho*p(2)+sigma_e2^2*P+N0);
%     
%     
%     h12_=[h1_ h2_];
%     %%%%caculation of p_c, f_c from (7),(8),(9),(10), with (15)
%     A=[h1_';h2_']*[h1_ h2_];%(10)
%     a_11=A(1,1);
%     a_12=A(1,2);
%     a_22=A(2,2);
%     
%     
%     
%     mu_1=(a_22-abs(a_12))/(a_11+a_22-2*abs(a_12));%(9)
%     mu_2=(a_11-abs(a_12))/(a_11+a_22-2*abs(a_12));
%     lam=(a_11*a_22-abs(a_12)^2)/(a_11+a_22-2*abs(a_12));%(8)
%     fc=(mu_1*h1_+mu_2*h2_*exp(-1i*angle(a_12)))/sqrt(lam);%*exp(1i*pi);
%     
%     norm(fc)^2;
% end






end


