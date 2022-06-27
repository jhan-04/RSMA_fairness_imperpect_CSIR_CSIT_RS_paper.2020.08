function [MA_p,tou_p, P1_p,P2_p, Pc_p,Rs_p]=RS_SDMA(P,h1,h2)
MA_p=1;

 rho=1-abs(h1'/norm(h1)*h2/norm(h2))^2;


Pc_p=0;
tou_p=1;

% mu=P/2+(1/(2*rho))*(1/norm(h2)^2+1/norm(h1)^2);
%  P1_p=mu-(1/rho)*(1/norm(h1)^2);
%  P2_p=mu-(1/rho)*(1/norm(h2)^2);
 a_max=100000000;
 a_min=-100000000;
 while (a_max-a_min)>0.0001
     a=(a_max+a_min)/2;
      P1_p=max(a-(1/rho)*(1/norm(h1)^2),0);
     P2_p=max(a-(1/rho)*(1/norm(h2)^2),0);
     
     if (P1_p+P2_p)>P
         a_max=a;
     else
         a_min=a;
     end
     
 end
 
 
Rs_p=log2(1+norm(h1)^2*rho*P1_p)+log2(1+norm(h2)^2*rho*P2_p);

end
    