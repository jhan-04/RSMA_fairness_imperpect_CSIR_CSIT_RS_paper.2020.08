
function [Rs,R1,R2,Rc]=cal_ach_R(P,t, P1,P2, Pc,fc,h1,h2,rho,beta)



H=[h1 h2];
w=H*inv(H'*H);
f1=w(:,1)/norm(w(:,1));
f2=w(:,2)/norm(w(:,2));

for sam=1:1:1000
    
    Nt=length(h1);
    e1=sqrt(beta)*(randn(Nt,1)+1i*randn(Nt,1))/sqrt(2);
    e2=sqrt(beta)*(randn(Nt,1)+1i*randn(Nt,1))/sqrt(2);
    hh1=h1+e1;
    hh2=h2+e2;
    
    
    Rc1=log2(1+abs(h1'*fc)^2*Pc/(1+abs(e1'*fc)^2*Pc+abs(hh1'*f1)^2*P1+abs(hh1'*f2)^2*P2));
    Rc2=log2(1+abs(h2'*fc)^2*Pc/(1+abs(e2'*fc)^2*Pc+abs(e2'*f1)^2*P1+abs(hh2'*f2)^2*P2)); 
    Rc_x(sam)=min(Rc1,Rc2);
    
    R1_x(sam)=log2(1+abs(h1'*f1)^2*P1/(1+abs(e1'*fc)^2*Pc+abs(e1'*f1)^2*P1+abs(hh1'*f2)^2*P2));
    R1_x2(sam)=log2(1+abs(h1'*f1)^2*P1/(1+abs(e1'*fc)^2*Pc+abs(e1'*f1)^2*P1+abs(e1'*f2)^2*P2));
    
    R2_x(sam)=log2(1+abs(h2'*f2)^2*P2/(1+abs(e2'*fc)^2*Pc+abs(hh2'*f1)^2*P1+abs(e2'*f2)^2*P2));
    R2_x2(sam)=log2(1+abs(h2'*f2)^2*P2/(1+abs(e2'*fc)^2*Pc+abs(e2'*f1)^2*P1+abs(e2'*f2)^2*P2));
    
    
    Rs_x(sam)=Rc_x(sam)+R1_x(sam)+R2_x(sam);
    
end

Rs=mean(Rs_x);
R1=mean(R1_x);
R2=mean(R2_x);
Rc=mean(Rc_x);


sum(abs(R1_x-R1_x2));
sum(abs(R2_x-R2_x2));


end