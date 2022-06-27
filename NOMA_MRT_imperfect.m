function [MA_x,t_x, P1_x,P2_x, Pc_x,Rs_x,fc_x]=NOMA_MRT_imperfect(P,h1,h2,rho,beta)

MA_x=2;
fc_x=h2/norm(h2);
N0=1+beta*P;
% Nt=2;
% P=10;
% ha=(randn(Nt,1)+1i*randn(Nt,1))/sqrt(2);%Rayleigh
% hb=(randn(Nt,1)+1i*randn(Nt,1))/sqrt(2);
% %     h1=ha;
% %     h2=ha;
% if norm(ha)>=norm(hb)
%     h1=ha;
%     h2=hb;
% else
%     h1=hb;
%     h2=ha;
% end
% rho=1-abs(h1'/norm(h1)*h2/norm(h2))^2;
% beta=0.8;
%


%%%case1  Rc1<Rc2
t_1=1;
R1_1=log2(N0+norm(h1)^2*P*t_1)-log2(N0);
Rc1_1=log2(1+norm(h1)^2*(1-rho)*P*(1-t_1)/(N0+norm(h1)^2*P*t_1));
Rc2_1=log2(1+norm(h2)^2*P*(1-t_1)/(N0+norm(h2)^2*(1-rho)*P*t_1));

Rs_1=Rc1_1+R1_1;


if Rc1_1>Rc2_1
    Rs_1=-10
end

%%%case2  Rc1==Rc2
t_2=(N0*(1-rho)*norm(h1)^2-norm(h2)^2)/((2*rho-rho^2)*norm(h1)^2*norm(h2)^2*P);
R1_2=log2(N0+norm(h1)^2*P*t_2)-log2(N0);
Rc1_2=log2(1+norm(h1)^2*(1-rho)*P*(1-t_2)/(N0+norm(h1)^2*P*t_2));
Rc2_2=log2(1+norm(h2)^2*P*(1-t_2)/(N0+norm(h2)^2*(1-rho)*P*t_2));
Rs_2=Rc1_2+R1_2;

if (0>t_2)||(t_2>1)
    Rs_2=-100;
end
    
    
    
    
    %%%case3  Rc1>Rc2
    
    a=-rho*norm(h1)^2*norm(h2)^2*P^2;
    b=P*(N0*(norm(h1)^2-rho*norm(h2)^2)+norm(h1)^2*norm(h2)^2*P);
    c=N0*(norm(h2)^2*P+N0);
    
    m=(1-rho)*norm(h2)^2*P;
    n=N0;
    syms t_4
    A=a*m;
    B=2*a*n;
    C=b*n-c*m;
    t3=double(solve((A*t_4^2+B*t_4+C)==0,t_4));
    
    tt31=(-B+sqrt(B^2-4*A*C))/(2*a);
    tt32=(-B-sqrt(B^2-4*A*C))/(2*a);
    
    if (0<=t3(1))&&(t3(1)<=1)
        t_3=double(t3(1));
    elseif (0<=t3(2))&&(t3(2)<=1)
        t_3=double(t3(2));
    else
        t_3=1;
    end
    Rs_3=log2((a*t_3^2+b*t_3+c)/(m*t_3+n))-log2(N0);
    
    R1_3=log2(N0+norm(h1)^2*P*t_3)-log2(N0);
    Rc1_3=log2(1+norm(h1)^2*(1-rho)*P*(1-t_3)/(N0+norm(h1)^2*P*t_3));
    Rc2_3=log2(1+norm(h2)^2*P*(1-t_3)/(N0+norm(h2)^2*(1-rho)*P*t_3));
    Rs_33=R1_3+Rc2_3;%Rs_33=Rs_3
    Rs_3-Rs_33;
    t_3;
    
    % syms t_4
    % Rs(t_4)=log2((a*t_4^2+b*t_4+c)/(m*t_4+n));
    % YY(t_4)=log2((N0+norm(h1)^2*P*t_4)*(1+norm(h2)^2*P*(1-t_4)/(N0+norm(h2)^2*(1-rho)*P*t_4)));
    % plot([0:0.1:1],2.^(YY([0:0.1:1])),'-r',[0:0.1:1],2.^(Rs([0:0.1:1])),':b');
    
    
    
    
    if Rc1_3<Rc2_3
        Rs_3=-10;
    end
    
    
    Rs=[Rs_1,Rs_2,Rs_3;Rc1_1,Rc2_1,0];
    t=[t_1,t_2,t_3,t3'];
    mama=max(Rs_1, Rs_2);
    ma=max(mama, Rs_3);
    if ma==Rs_1
        Rs_x=Rs_1;
        t_x=t_1;
    end
    if ma==Rs_2
        Rs_x=Rs_2;
        t_x=t_2;
    end
    if ma==Rs_3
        Rs_x=Rs_3;
        t_x=t_1;
    end
    
    
    P2_x=0;
    Pc_x=P*(1-t_x);
    P1_x=P*t_x;
    
    
    
    
    
end
