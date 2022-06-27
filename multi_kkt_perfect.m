function [MA_x,t_x, P1_x,P2_x, Pc_x,Rs_x,fc_x]=multi_kkt_perfect(P,h1,h2,rho)


MA_x=4;
P2_x=0;
P1_x=0;
Pc_x=P;
t_x=0;
Nt=length(h1);


if rho==0
    
    fc_x=h1/norm(h1);
else
    h1_=h1;
    h2_=h2;
    
    
    h12_=[h1_ h2_];
    %%%%caculation of p_c, f_c from (7),(8),(9),(10), with (15)
    A=[h1_';h2_']*[h1_ h2_];%(10)
    a_11=A(1,1);
    a_12=A(1,2);
    a_22=A(2,2);
    
    
    
    mu_1=(a_22-abs(a_12))/(a_11+a_22-2*abs(a_12));%(9)
    mu_2=(a_11-abs(a_12))/(a_11+a_22-2*abs(a_12));
    lam=(a_11*a_22-abs(a_12)^2)/(a_11+a_22-2*abs(a_12));%(8)
    fc_x=(mu_1*h1_+mu_2*h2_*exp(-1i*angle(a_12)))/sqrt(lam);%*exp(1i*pi);
    
    norm(fc_x)^2;
end

Rs_x=min(log2(1+abs(h1'*fc_x)^2*P),log2(1+abs(h2'*fc_x)^2*P));




end