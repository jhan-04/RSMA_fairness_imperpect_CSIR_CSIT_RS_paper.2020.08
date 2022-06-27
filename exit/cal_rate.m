function [Rs]=cal_rate(P,tou_x, P1_x,P2_x, Pc_x,fc_x,f1_x,f2_x,hh1,hh2, beta)


%hh1=h1+e1
%hh2=h2+e2
    Rc1=log2(1+abs(hh1'*fc_x)^2*Pc_x/(1+abs(hh1'*f1_x)^2*P1_x+abs(hh1'*f2_x)^2*P2_x));
  
    Rc2=log2(1+abs(hh2'*fc_x)^2*Pc_x/(1+abs(hh2'*f1_x)^2*P1_x+abs(hh2'*f2_x)^2*P2_x));
    
    Rc_x=min(Rc1,Rc2);
    R1_x=log2(1+abs(hh1'*f1_x)^2*P1_x/(1+abs(hh1'*f2_x)^2*P2_x));
   
    R2_x=log2(1+abs(hh2'*f2_x)^2*P2_x/(1+abs(hh2'*f1_x)^2*P1_x));
    Rs=Rc_x+R1_x+R2_x;


end