%Rate-Splitting Unifying SDMA, OMA, NOMA, and Multicasting in MISO Broadcast Channel_ A Simple Two-User Rate Analysis
%Fig 3 (a)
function [MA_x,tou_x, P1_x,P2_x, Pc_x,Rs_x,fc_x]=NOMA_kkt_imperfect(P,h1,h2,rho,beta)


% rho=1-abs(h1'/norm(h1)*h2/norm(h2))^2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%case2)
v=0;
for t=0:1/500:1
    
    v=v+1;
    h1_=h1/sqrt((1+beta*P)+norm(h1)^2*rho*P*t);
    h2_=h2/sqrt(1+beta*P);
    A=[h1_';h2_']*[h1_ h2_];%(10)
    
    if rho==0
        
        fc=h1/norm(h1);
    else
        
        a_11=A(1,1);
        a_12=A(1,2);
        a_22=A(2,2);
        
        mu_1=(a_22-abs(a_12))/(a_11+a_22-2*abs(a_12));%(9)
        mu_2=(a_11-abs(a_12))/(a_11+a_22-2*abs(a_12));
        lam=(a_11*a_22-abs(a_12)^2)/(a_11+a_22-2*abs(a_12));%(8)
        fc=(mu_1*h1_+mu_2*h2_*exp(-1i*angle(a_12)))/sqrt(lam);
        
    end
    
    b1=P*abs(h2_'*fc)^2;
    a1=norm(h1)^2*rho*P;
    %a1=norm(h1)^2*P;
    R(v)=log2((1+beta*P)+a1*t)+log2(1+b1*(1- t))-log2(1+beta*P);%R1+Rc(Rc1==Rc2)
    Rc1(v)=log2(1+abs(h1_'*fc)^2*(1- t)*P);
    Rc2(v)=log2(1+abs(h2_'*fc)^2*(1- t)*P);
    
    fc_22(:,v)=fc;
end

t=0:1/500:1;
%plot(t,R)
Rs_2=max(R);%0<=t<1
k=find(Rs_2==R);
t2=t(k);
fc_2=fc_22(:,k);


    
    tou_x=t2;
    Rs_x=Rs_2;
    P2_x=0;
    P1_x=P*t2;
    Pc_x=(1-t2)*P;
    fc_x=fc_2;
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if tou_x==1
    MA_x=1;%SDMA
end
if (P1_x~=0)&&(P2_x==0)&&(Pc_x~=0)
    MA_x=2;%NOMA
end
if (P1_x~=0)&&(P2_x==0)&&(Pc_x==0)
    MA_x=3;%OMA
end
if (P2_x==0)&&(P1_x==0)%tou(i,j)==0
    MA_x=4;%multicasting
end
if (P1_x~=0)&&(P2_x~=0)&&(Pc_x~=0)
    MA_x=5;%RS
end




end



