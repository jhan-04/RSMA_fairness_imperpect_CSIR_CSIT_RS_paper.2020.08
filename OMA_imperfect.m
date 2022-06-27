function [MA_x,t_x, P1_x,P2_x, Pc_x,Rs_x,fc_x]=OMA_imperfect(P,h1,h2,rho,beta)


MA_x=3;

P2_x=0;
Pc_x=0;
P1_x=P;

t_x=1;
Nt=length(h1);
fc_x=zeros(Nt,1);

Rs_x=log2(1+beta*P+norm(h1)^2*rho*P1_x)-log2(1+beta*P);

end 
    