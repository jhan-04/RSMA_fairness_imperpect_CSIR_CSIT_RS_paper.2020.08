
function  [MA_x,tou_x, P1_x,P2_x, Pc_x,Rs_x,Rs_noma_mul_oma]=RS_summaxrate(P,h1,h2)

 rho=1-abs(h1'/norm(h1)*h2/norm(h2))^2;



h1=1/sqrt(2)*[1;1];
h2=(gam)/sqrt(2)*[1;exp(-1i*theta)];
h1_=h1/norm(h1);
h2_=h2/norm(h2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A=h1_*h1_';
B=h2_*h2_';
%%%%%%%%%%%%%%%%%%%
cvx_begin quiet
variable F(2,2) hermitian
variable x(1)
minimize(-x)
subject to
trace(A*F)>=x;
trace(B*F)>=x;
trace(F)==1;
F== semidefinite(2);
cvx_end
%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%
F;




[U,W,Z] = svds(F);%V=U*W*Z'
obj_max=-100;
n=length(W);
for m=1:10000
    r=sqrt(1/2)*(randn(2,1)+1i*randn(2,1));%sqrt(var/2)*(randn(1,N)+1i*randn(1,N))
    v_=U*sqrt(W)*r;
    v_=v_/sqrt(v_'*v_);
    obj=min(v_'*A*v_,v_'*B*v_);
    if obj>=obj_max
        obj_max=obj;
        v_max=v_;
    end
end
f_c=v_max;
(f_c'*A*f_c)-(f_c'*B*f_c);

%%%%%%%P1 ~=0,P2 ~=0

if (f_c'*A*f_c)>=(f_c'*B*f_c)
    ha=h1;
    hb=h2;
    Gam=(1/rho)*(1/norm(hb)^2-1/norm(ha)^2);
    b=norm(ha)^2*rho*P/2;
    a=1+(Gam/P)*b;
    d=norm(hb)^2*rho*P/2-abs(hb'*f_c)^2*P;
    c=1-(Gam/P)*d+abs(hb'*f_c)^2*(P-Gam);
    
else
    ha=h2;
    hb=h1;
    Gam=(1/rho)*(1/norm(hb)^2-1/norm(ha)^2);
    b=norm(ha)^2*rho*P/2;
    a=1+(Gam/P)*b;
    d=norm(hb)^2*rho*P/2-abs(hb'*f_c)^2*P;
    c=1-(Gam/P)*d+abs(hb'*f_c)^2*(P-Gam);
end
if (b*d)>0
    
    t12=1;
    
    Rs_121=log2(a*c+(a*d+b*c)*t12+b*d*t12^2);
    t12=0;
    Rs_120=log2(a*c+(a*d+b*c)*t12+b*d*t12^2);
    if  Rs_121>Rs_120
        t12=1;
       Rs_12=log2(a*c+(a*d+b*c)*t12+b*d*t12^2);
    else
         t12=0;
       Rs_12=log2(a*c+(a*d+b*c)*t12+b*d*t12^2);
    end
    
    
elseif (b*d)==0
    t12=((a*d+b*c)+abs(a*d+b*c))/(2*a*d+2*b*c);
    
    Rs_12=log2(a*c+(a*d+b*c)*t12+b*d*t12^2);
    1111
else
    
    t12=min(-a/(2*b)-c/(2*d),1);
    
    Rs_12=log2(a*c+(a*d+b*c)*t12+b*d*t12^2);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%1


if (t12*P >= (1/rho)*(1/norm(h2)^2-1/norm(h1)^2)) %%P2~=0 P2가 0이 아니면 당연히 P1도 0이 아니다.
    Rs_1=Rs_12;
    tou_1=t12;
    P1_1=t12*P/2+(-1/norm(h1)^2+1/norm(h2)^2)/(2*rho);
    P2_1=t12*P/2+(1/norm(h1)^2-1/norm(h2)^2)/(2*rho);
    Pc_1=(1-t12)*P;
    % %     Rs_x=Rs_12;
    % %     tou_x=t12;
    % %     P1_x=t12*P/2+(-1/norm(h1)^2+1/norm(h2)^2)/(2*rho);
    % %     P2_x=t12*P/2+(1/norm(h1)^2-1/norm(h2)^2)/(2*rho);
    % %     Pc_x=(1-t12)*P;
    
    
else
    Rs_1=-100;
end
v=0;


for t=0:1/200:1
    v=v+1;
    h1_=h1/sqrt(1+norm(h1)^2*rho*P*t);
    % h1_=h1/sqrt(1+norm(h1)^2*P*t);
    h2_=h2;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    A=h1_*h1_';
    B=h2_*h2_';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cvx_begin quiet
    variable F(2,2) hermitian
    variable x(1)
    minimize(-x)
    subject to
    trace(A*F)>=x;
    trace(B*F)>=x;
    trace(F)==1;
    F== semidefinite(2);
    cvx_end
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    trace(A*F);
    trace(B*F);
    F;
    
    [U,W,Z] = svds(F);%V=U*W*Z'
    obj_max=-100;
    n=length(W);
    for m=1:10000
        r=sqrt(1/2)*(randn(n,1)+1i*randn(n,1));%sqrt(var/2)*(randn(1,N)+1i*randn(1,N))
        v_=U*sqrt(W)*r;
        v_=v_/sqrt(v_'*v_);
        obj=min(v_'*A*v_,v_'*B*v_);
        if obj>=obj_max
            obj_max=obj;
            v_max=v_;
        end
    end
    f_c=v_max;
    norm(f_c)^2;
    (f_c'*A*f_c)-(f_c'*B*f_c);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Rc1(v)=log2(1+abs(h1_'*f_c)^2*(1- t)*P);
    Rc2(v)=log2(1+abs(h2_'*f_c)^2*(1- t)*P);
    
    Rc=min(log2(1+abs(h1_'*f_c)^2*(1- t)*P),log2(1+abs(h2_'*f_c)^2*(1- t)*P));
    log2(1+abs(h1_'*f_c)*(1- t)*P);
    log2(1+abs(h2_'*f_c)*(1- t)*P);
    R(v)=  Rc+ log2(1+norm(h1)^2*rho*P*t);
    %         b1=P*abs(h2'*f_c)^2;
    %         a1=norm(h1)^2*P;
    %         R(v)=log2(1+a1*t)+log2(1+b1*(1- t));
    
    abs(h1'*f_c)^2;
    abs(h2'*f_c)^2;
end
R;
t=0:1/200:1;
Rs_noma_mul_oma=max(R);%0<=t<1
k=find(Rs_noma_mul_oma==R);
t1=t(k);
RC1_m=Rc1(k);
RC2_m=Rc2(k);



maax= max(Rs_1,Rs_noma_mul_oma);

if maax==Rs_noma_mul_oma
    tou_x=t1;
    Rs_x=Rs_noma_mul_oma;
    P2_x=0;
    P1_x=P*t1;
    Pc_x=(1-t1)*P;
else
    Rs_x=Rs_12;
    tou_x=t12;
    P1_x=t12*P/2+(-1/norm(h1)^2+1/norm(h2)^2)/(2*rho);
    P2_x=t12*P/2+(1/norm(h1)^2-1/norm(h2)^2)/(2*rho);
    Pc_x=(1-t12)*P;
    
end


%end


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







