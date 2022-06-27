 %Rate-Splitting Unifying SDMA, OMA, NOMA, and Multicasting in MISO Broadcast Channel_ A Simple Two-User Rate Analysis
%Fig 3 (a)
clear all
clc
nt=2;%number of transmitter antenna
P=100;
%tilde is ~
%gam=rand();%channel strength diversity
%rho=1-(abs(h1_'*h2_))^2;%angle of between user channel direction
i=0;
for gam_dB=-20:0.1:0
    
    gam=10^(gam_dB/20);
    i=i+1;
    j=0;
    k=0;
    for  rho=0.01:0.1:1
    % for rho=0:0
        j=j+1;
        
        
        theta=acos(1-2*rho);
        Theta(i,j)=acos(1-2*rho);
        
        
        h1=1/sqrt(2)*[1;1];
        h2=(gam)/sqrt(2)*[1;exp(-1i*theta)];
        h1_=h1;
        h2_=h2/gam;
        
        
        h12_=[h1_ h2_];
        %%%%caculation of p_c, f_c from (7),(8),(9),(10), with (15)
        A=[h1_';h2_']*[h1_ h2_];%(10)
        a_11=A(1,1);
        a_12=A(1,2);
        a_22=A(2,2);
        
        
        
        mu_1=(a_22-abs(a_12))/(a_11+a_22-2*abs(a_12));%(9)
        mu_2=(a_11-abs(a_12))/(a_11+a_22-2*abs(a_12));
        lam=(a_11*a_22-abs(a_12)^2)/(a_11+a_22-2*abs(a_12));%(8)
        f_c=(mu_1*h1_+mu_2*h2_*exp(-1i*angle(a_12)))/sqrt(lam);%*exp(1i*pi);
        
        
        
        
        %%%%%%caculation of
        
        %%%%%%%P1 ~=0,P2 ~=0
        
        Gam=(1/rho)*(1/norm(h2)^2-1/norm(h1)^2);
        b=norm(h1)^2*rho*P/2;
        a=1+(Gam/P)*b;
        d=norm(h2)^2*rho*P/2-abs(h2'*f_c)^2*P;
        c=1-(Gam/P)*d+abs(h2'*f_c)^2*(P-Gam);
        
        
        
        
        t12(i,j)=min(-a/(2*b)-c/(2*d),1);
        
        Rs_12(i,j)=log2(a*c+(a*d+b*c)*t12(i,j)+b*d*t12(i,j)^2);
        
        
      
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%1
        
        
        if (t12(i,j)*P >= Gam) %%P2~=0 P2가 0이 아니면 당연히 P1도 0이 아니다.
            Rs_x(i,j)=Rs_12(i,j);
            tou_x(i,j)=t12(i,j);
            P1_x=t12(i,j)*P/2+(-1/norm(h1)^2+1/norm(h2)^2)/(2*rho);
            P2_x=t12(i,j)*P/2+(1/norm(h1)^2-1/norm(h2)^2)/(2*rho);
            Pc_x=(1-t12(i,j))*P;
            
            % dis(i,j)=abs(h1'*f_c)^2/(1+norm(h1)^2*rho*P1_x)-abs(h2'*f_c)^2/(1+norm(h2)^2*rho*P2_x);
        else
            v=0;
            for t=0:1/1000:1
                v=v+1;
                %h1_=h1/sqrt(1+norm(h1)^2*rho*P*t);
                
                h1_=h1/sqrt(1+norm(h1)^2*P*t);
                h2_=h2;
                A=[h1_';h2_']*[h1_ h2_];%(10)
                a_11=A(1,1);
                a_12=A(1,2);
                a_22=A(2,2);
                
                mu_1=(a_22-abs(a_12))/(a_11+a_22-2*abs(a_12));%(9)
                mu_2=(a_11-abs(a_12))/(a_11+a_22-2*abs(a_12));
                lam=(a_11*a_22-abs(a_12)^2)/(a_11+a_22-2*abs(a_12));%(8)
                f_c=(mu_1*h1_+mu_2*h2_*exp(-1i*angle(a_12)))/sqrt(lam);
                    
                 
                b1=P*abs(h2'*f_c)^2;
                %a1=norm(h1)^2*rho*P;
                a1=norm(h1)^2*P;
                R(v)=log2(1+a1*t)+log2(1+b1*(1- t));
               %R(v)=log2(1+norm(h1)^2*rho*t*P)+log2(1+abs(h2'*f_c));
                
            end

% syms t
% h1_=h1/sqrt(1+norm(h1)^2*rho*P*t);
% h2_=h2;
% 
% equ=diff(log2(1+norm(h1)^2*rho*P*t)+log2(1+P*abs(h2'*{(a_22-abs(a_12)/sqrt(1+norm(h1)^2*rho*P*t))*h1/sqrt(1+norm(h1)^2*rho*P*t)})^2),t);



            
            t=0:1/1000:1;
            Rs_1(i,j)=max(R);
            k=find(Rs_1(i,j)==R);
            t1(i,j)=t(k);
            
            tou_x(i,j)=t1(i,j);
            P2_x=0;
            P1_x=P*t1(i,j);
            Pc_x=(1-t1(i,j))*P;
            Rs_x(i,j)=Rs_1(i,j);
            
            
            
%             if tou_x(i,j)*P/2+(-1/norm(h1)^2+1/norm(h2)^2)/(2*rho) <0
%                 11
%                 
%                 Rs_c(i,j)=min(log2(1+P*abs(h1'*f_c)^2),log2(1+P*abs(h2'*f_c)^2));
%                 tou_x(i,j)=0;
%                 P2_x=0;
%                 P1_x=0;
%                 Pc_x=P;
%                 Rs_x(i,j)=Rs_c(i,j);
%                 
%             end
            
            %dis(i,j)=abs(h1'*f_c)^2/(norm(h1)^2)-abs(h2'*f_c)^2/(norm(h2)^2);
            %dis(i,j)=abs(h1'*f_c)^2/(1+norm(h1)^2*rho*P1_x)-abs(h2'*f_c)^2/(1+norm(h2)^2*rho*P2_x);
            
        end
        
        
        if tou_x(i,j)==1
            MA_x(i,j)=1;%SDMA
        end
        if (P1_x~=0)&&(P2_x==0)&&(Pc_x~=0)
            MA_x(i,j)=2;%NOMA
        end
        if (P1_x~=0)&&(P2_x==0)&&(Pc_x==0)
            MA_x(i,j)=3;%OMA
        end
        if (P2_x==0)&&(P1_x==0)%tou(i,j)==0
            MA_x(i,j)=4;%multicasting
        end
        if (P1_x~=0)&&(P2_x~=0)&&(Pc_x~=0)
            MA_x(i,j)=5;%RS
        end
        
        
        
        
    end
    
end



figure(1)
x=0.01:0.1:1
y=-20:0.1:0;
[X Y]=meshgrid(x,y);
% contourf(X,Y,tou,'ShowText','on')


%contourf(X,Y,MA_x)
contourf(X,Y,tou_x)
contourcbar

xlabel('rho')
ylabel('channel strength disparity  [dB]')
figure(2)
contourf(X,Y,MA_y,5)
%
