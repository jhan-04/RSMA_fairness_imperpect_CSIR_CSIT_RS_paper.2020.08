function [Rs_x,Rs_x1,Rs_x2]=OMA_MRT_imperfect(P,h1,h2,rho,beta)

MA_x=3;

Nt=length(h1);fc_x=zeros(Nt,1);

P2_x=0;
Pc_x=0;
P1_x=P;

t_x=1;

Rs_x1=1/2*log2(1+norm(h1)^2*P/(1+beta*P));
Rs_x2=1/2*log2(1+norm(h2)^2*P/(1+beta*P));
Rs_x=1/2*log2(1+norm(h1)^2*P/(1+beta*P))+1/2*log2(1+norm(h2)^2*P/(1+beta*P));

end
    