
function  [fc]=Find_fc_kkt(P,p,h1,h2,rho,N0, sigma_e1,sigma_e2)




if rho==0
    
    fc=h1/norm(h1);
else
h1_=h1/sqrt(norm(h1)^2*rho*p(1)+sigma_e1^2*P+N0);
h2_=h2/sqrt(norm(h2)^2*rho*p(2)+sigma_e2^2*P+N0);
    
    
    h12_=[h1_ h2_];
    %%%%caculation of p_c, f_c from (7),(8),(9),(10), with (15)
    A=[h1_';h2_']*[h1_ h2_];%(10)
    a_11=A(1,1);
    a_12=A(1,2);
    a_22=A(2,2);
    
    
    
    mu_1=(a_22-abs(a_12))/(a_11+a_22-2*abs(a_12));%(9)
    mu_2=(a_11-abs(a_12))/(a_11+a_22-2*abs(a_12));
    lam=(a_11*a_22-abs(a_12)^2)/(a_11+a_22-2*abs(a_12));%(8)
    fc=(mu_1*h1_+mu_2*h2_*exp(-1i*angle(a_12)))/sqrt(lam);%*exp(1i*pi);
    
    norm(fc)^2;
end






end


