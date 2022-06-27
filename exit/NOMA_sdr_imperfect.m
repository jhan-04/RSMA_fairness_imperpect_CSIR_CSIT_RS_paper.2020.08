
function  [MA_x,t_x, P1_x,P2_x, Pc_x,Rs_x,fc_x]=NOMA_sdr_imperfect(P,h1,h2,rho,beta)



Nt=length(h1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%CASE(1)P1>0,P2>0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h1_=h1/norm(h1);
h2_=h2/norm(h2);
A=h1_*h1_';
B=h2_*h2_';




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%P1>=0,P2=0,(P1=tP,P2=0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%CASE(2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%CASE(2)-1,Rc1>Rc2


fc_2_1=h2/norm(h2);
PP=P/(1+beta*P);
e=-PP^2*norm(h1)^2*norm(h2)^2*rho;
f=PP*(norm(h1)^2*rho-norm(h2)^2+norm(h1)^2*rho*norm(h2)^2*PP);
g=1+norm(h2)^2*PP;
if e<0
    t2_1=max(0,min(-f/(2*e),1));
else
    t2_1=0;
end

Rs_2_1=log2(e*t2_1^2+f*t2_1+g);
%%가정이 틀릴경우
if log2(1+norm(h1)^2*(1-rho)*(1-t2_1)*P/((1+beta*P)+P*norm(h1)^2*rho*t2_1))<=log2(1+norm(h2)^2*(1- t2_1)*P/(1+beta*P))
    Rs_2_1=-100;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%CASE(2)-3 ,Rc1==Rc2
U=100000;
L=1;
T=(U+L)/2;

k_min=(1+beta*P);
k_max=(1+beta*P)+norm(h1)^2*rho*P;
k_region=[k_min:(k_max- k_min)/100:k_max];
syms k real
QQ(k)=k+k*(norm(h2)^2*(norm(h1)^2*P*rho+(1+beta*P)-k)/(k*norm(h2)^2+(1+beta*P)*norm(h1)^2-2*sqrt(k*(1+beta*P))*norm(h1)*norm(h2)*sqrt(1-rho)));

while abs((U-L))>0.01
    
    T=(U+L)/2;
    
    cond1=QQ(k)>=T;
    cond2=(1+beta*P)<=k;
    cond3=k<=(1+beta*P)+norm(h1)^2*rho*P;
    conds=[cond1 cond2 cond3];
    sol = solve(conds, k, 'Real', true,'ReturnConditions', true, 'IgnoreAnalyticConstraints', true);
    sol.k;
    sol.parameters;
    sol.conditions;
    
    if isempty(sol.conditions)
        U=T;
    else
        condWithValues = subs(sol.conditions(1), sol.k(1),k_region );
        feasible=isAlways(condWithValues);
        if sum(feasible==1)==0
            U=T;
        else
            L=T;
            d=find(feasible==1);
            k_min=k_region(min(d));
            k_max=k_region(max(d));
            k_region=[k_min:(k_max-k_min)/100:k_max];
            if k_min==k_max
                K=k_max;
            end
        end
        
    end
    
    if k_min==k_max
        break;
    end
    K= mean(k_region);
    
end

t2_2=round((K-(1+beta*P))/(norm(h1)^2*rho*P),4);
%%%%%
h1_=h1/sqrt((1+beta*P)+norm(h1)^2*rho*P*t2_2);
h2_=h2/sqrt(1+beta*P);
A=[h1_';h2_']*[h1_ h2_];%(10)
a_11=A(1,1);
a_12=A(1,2);
a_22=A(2,2);
mu_1=(a_22-abs(a_12))/(a_11+a_22-2*abs(a_12));%(9)
mu_2=(a_11-abs(a_12))/(a_11+a_22-2*abs(a_12));
lam=(a_11*a_22-abs(a_12)^2)/(a_11+a_22-2*abs(a_12));%(8)
fc_2_2=(mu_1*h1_+mu_2*h2_*exp(-1i*angle(a_12)))/sqrt(lam);
%%%%%%%
Rs_2_2=log2(double(QQ(K)))-log2(1+beta*P);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%CASE(2)-3가지 경우 중 좋은 것

a=max(Rs_2_1,Rs_2_2);
if a==Rs_2_1
    t2=t2_1;
    Rs_2=Rs_2_1;
    fc_2=fc_2_1;
    
end
if a==Rs_2_2
    t2=t2_2;
    Rs_2=Rs_2_2;
    fc_2=fc_2_2;
    
end


    t_x=t2;
    Rs_x=Rs_2;
    P2_x=0;
    P1_x=P*t2;
    Pc_x=(1-t2)*P;
    fc_x=fc_2;
 


if t_x==1
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


