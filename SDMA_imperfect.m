function [MA_x,t_x, P1_x,P2_x, Pc_x,Rs_x,fc_x]=SDMA_imperfect(P,h1,h2,rho,beta)


MA_x=1;
Pc_x=0;
t_x=1;
Nt=length(h1);
fc_x=0*rand(Nt,1);
% mu=P/2+(1/(2*rho))*(1/norm(h2)^2+1/norm(h1)^2);
%  P1_p=mu-(1/rho)*(1/norm(h1)^2);
%  P2_p=mu-(1/rho)*(1/norm(h2)^2);
 a_max=100000000;
 a_min=-100000000;
 while (a_max-a_min)>0.0001
     a=(a_max+a_min)/2;
      P1_x=max(a-((1+beta*P)/rho)*(1/norm(h1)^2),0);
     P2_x=max(a-((1+beta*P)/rho)*(1/norm(h2)^2),0);
     
     if (P1_x+P2_x)>P
         a_max=a;
     else
         a_min=a;
     end
     
 end
 
 
Rs_x=log2(1+norm(h1)^2*rho*P1_x/(1+beta*P))+log2(1+norm(h2)^2*rho*P2_x/(1+beta*P));

    
end