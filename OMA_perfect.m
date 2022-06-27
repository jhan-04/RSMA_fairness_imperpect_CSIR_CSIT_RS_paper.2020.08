function [MA_x,t_x, P1_x,P2_x, Pc_x,Rs_x,fc_x]=OMA_perfect(P,h1,h2,rho)


MA_x=3;

P2_x=0;
Pc_x=0;
P1_x=P;

t_x=1;

fc_x=zeros(2,1);
Rs_x=log2(1+norm(h1)^2*P*rho);

end
    