
function  [MA_x,t_x, P1_x,P2_x, Pc_x,Rs_x,fc_x]=NOMA_sdr_perfect(P,h1,h2,rho)

if rho==0
    
    fc_x=h1/norm(h1);
    
    t_x=0;
    h1_=h1/sqrt(1+norm(h1)^2*rho*P*t2);
    
    h2_=h2;
    Rc1=log2(1+abs(h1_'*f_c)^2*(1- t)*P);
    Rc2=log2(1+abs(h2_'*f_c)^2*(1- t)*P);
    Rs_x=min(Rc1,Rc2);
    
    P1_x=0;
    P2_x=0;
    Pc_x=P;
    
else
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%P1>=0,P2=0,(P1=tP,P2=0)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%CASE(2)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%CASE(2)-1,Rc1>Rc2
    
    
    fc_2_1=h2/norm(h2);
    e=-P^2*norm(h1)^2*norm(h2)^2*rho;
    f=P*(norm(h1)^2*rho-norm(h2)^2+norm(h1)^2*rho*norm(h2)^2*P);
    g=1+norm(h2)^2*P;
    if e<0
        t_2_1=max(0,min(-f/(2*e),1));
    else
        t_2_1=0;
    end
    
    Rs_2_1=log2(e*t_2_1^2+f*t_2_1+g);
    %%가정이 틀릴경우
    if log2(1+norm(h1)^2*(1-rho)*(1-t_2_1)*P/(1+P*norm(h1)^2*rho*t_2_1))<=log2(1+norm(h2)^2*(1- t_2_1)*P)
        Rs_2_1=-10;
    end
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%CASE(2)-3 ,Rc1==Rc2
    U=1000;
    L=1;
    T=(U+L)/2;
  
    k_min=1;
    k_max=1+norm(h1)^2*rho*P;
    k_region=[ k_min:(k_max- k_min)/100:k_max];
    syms k real
    QQ(k)=k+k*(norm(h2)^2*(norm(h1)^2*P*rho+1-k)/(k*norm(h2)^2+norm(h1)^2-2*sqrt(k)*norm(h1)*norm(h2)*sqrt(1-rho)));
    
    while abs((U-L)/L)>0.01
        
        T=(U+L)/2;
        
        cond1=QQ(k)>=T;
        cond2=1<=k;
        cond3=k<=1+norm(h1)^2*rho*P;
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
    
    t2_2=(K-1)/(norm(h1)^2*rho*P);  
    %%%%%
    h1_=h1/sqrt(1+norm(h1)^2*rho*P*t2_2);
    h2_=h2;
    A=[h1_';h2_']*[h1_ h2_];%(10) 
    a_11=A(1,1);
    a_12=A(1,2);
    a_22=A(2,2);  
    mu_1=(a_22-abs(a_12))/(a_11+a_22-2*abs(a_12));%(9)
    mu_2=(a_11-abs(a_12))/(a_11+a_22-2*abs(a_12));
    lam=(a_11*a_22-abs(a_12)^2)/(a_11+a_22-2*abs(a_12));%(8)
    fc_2_2=(mu_1*h1_+mu_2*h2_*exp(-1i*angle(a_12)))/sqrt(lam);
    %%%%%%%
    Rs_2_2=log2(double(QQ(K)));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%CASE(2)-3가지 경우 중 좋은 것
    
    a=max(Rs_2_1,Rs_2_2);
    if a==Rs_2_1
        t_2=t_2_1;
        Rs_2=Rs_2_1;
        fc_2=fc_2_1;
    end
    if a==Rs_2_2
        t_2=t2_2;
        Rs_2=Rs_2_2;
        fc_2=fc_2_2;
    end
    

        t_x=t_2;
        Rs_x=Rs_2;
        P2_x=0;
        P1_x=P*t_2;
        Pc_x=(1-t_2)*P;
        fc_x=fc_2;

end%rho==0

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


