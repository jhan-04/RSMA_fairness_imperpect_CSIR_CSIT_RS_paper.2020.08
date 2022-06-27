function [MA_x,t_x, P1_x,P2_x, Pc_x,Rs_x,fc_x]=OMA_MRT_perfect(P,h1,h2,rho)

MA_x=3;

Nt=length(h1);fc_x=zeros(Nt,1);

P2_x=0;
Pc_x=0;
P1_x=P;

t_x=1;


Rs_x=1/2*log2(1+norm(h1)^2*P)+1/2*log2(1+norm(h2)^2*P);

end
    