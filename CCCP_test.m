clc
clear all

Nt=2;%nuber of transmitter antenna
K=2;%number of user


SNR=[10];
error=sqrt(0.05);
sigma_e=error+zeros(K,1);
beta=(error)^2;
N0=1;
%i=1;
j=1;
num_interation=10;

for i=1:1
    i
    

H_h1=(randn(Nt,K)+1i*randn(Nt,K))/sqrt(2);%Rayleigh h1
H_h=sqrt(1-beta)*H_h1;
if norm(H_h(:,1))>=norm(H_h(:,2))
    h1(:,i)=H_h(:,1); h2(:,i)=H_h(:,2);
else
    h1(:,i)=H_h(:,2); h2(:,i)=H_h(:,1);
end
H_h(:,1)=h1(:,i);
H_h(:,2)=h2(:,i);


snr=SNR(1);
Pt=10^(snr/10);

[A,B,C,D]=cal_ABCD(Nt,K,Pt,N0,H_h,sigma_e);


[p_opt,result_cccp(i,:),GMI_cccp(i,:)]=GMI_only_cccp(Nt,K,H_h,Pt,N0,sigma_e,num_interation);
[p_opt2,result_cccp2(i,:),GMI_cccp2(i,:)]=GMI_only_cccp_ver2(Nt,K,H_h,Pt,N0,sigma_e,num_interation);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%CCCP and SDR %%%%%%%%%%%%%%

%initialization
b0(:,1)=(rand(K,1));%randn(M,1)*0+5;%
d0(:,1)=(rand(K,1));%randn(M,1)*0+5;%
ss=1;
result_SDR(1)=0;
GMI_SDR(1)=0;
stop=0;
%repeat
while ~stop
    
    ss=ss+1;
    if b0(1,ss-1)~=b0(1,ss-1)
        b0(:,ss-1)=(rand(K,1));
        d0(:,ss-1)=(rand(K,1));
    end
    
    cvx_begin quiet
    %cvx_begin
    variable X(Nt*(K+1),Nt*(K+1)) complex semidefinite
    variable Lc(1)
    variable a(K,1)
    variable b(K,1)
    variable c(K,1)
    variable d(K,1)
    xx=Lc+sum(c)-sum(d);
    maximize(xx)
    subject to
    for k=1:K
        a(k)-b(k)>=Lc;
        trace(A(:,:,k)*X)>=exp(a(k));
        trace(C(:,:,k)*X)>=exp(c(k));
        trace(B(:,:,k)*X)<=exp(b0(k,ss-1))*(b(k)-b0(k,ss-1)+1);
        trace(D(:,:,k)*X)<=exp(d0(k,ss-1))*(d(k)-d0(k,ss-1)+1);
    end
    %trace(X)<=Pt;
    X== hermitian_semidefinite(Nt*(K+1));
    
    cvx_end
    
    b0(:,ss)=b;
    d0(:,ss)=d;
    result_SDR(i,ss)=Lc+sum(c)-sum(d);
    
    
    if result_SDR(i,ss-1)>result_SDR(i,ss)&&ss>2
        fprintf('NOT monotonically increasing %d in Rs \n',result_SDR(i,ss-1)-result_SDR(i,ss))
        result_SDR
    end
    %
    [p_max]=find_p(K,Pt,A,B,C,D,X);%norm(p_max(:,i,j))^2;
        [GMI_SDR(i,ss),GMI_c,GMI_p]=cal_GMI(K,A,B,C,D,p_max);
        [GMI_SDR2(i,ss)]=cal_GMI_withX(K,A,B,C,D,X);
    if ss>num_interation%abs(GMI_SDR(ss)-GMI_SDR(ss-1))<ee
        stop=1;
    end
    
end

fprintf('iteration with SDR and cccp  %d in Rs \n',ss-2)
% %GMI==result/log(2)
result_SDR(i,ss);


end
%%%%%%%%%%%%%%%

figure
hold on
plot([0:num_interation],mean(result_SDR),'r-*')
plot([0:num_interation],mean(result_cccp),'b--*')
plot([0:num_interation],mean(result_cccp2),'g--*')
legend('SDR and CCCP','only CCCP')
title('lower bound')

figure
hold on
plot([0:num_interation],mean(GMI_SDR),'r-*')
plot([0:num_interation],mean(GMI_cccp),'b--*')
plot([0:num_interation],mean(GMI_cccp2),'g--*')
legend('SDR and CCCP','only CCCP')
title('GMI')