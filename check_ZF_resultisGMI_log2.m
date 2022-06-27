clear all
clc

%%%%ZF_확인
Nt=2;%nuber of transmitter antenna
M=2;%number of user
SNR=20;%[5:5:35];
%Pt=10^(SNR/10);%transmit signal power
N0=1;%gaussina noise

sigma_e1=0.2; sigma_e2=0.2; sigma_e3=0.2;
sigma_e=randn(M,1)*0+sigma_e1;%M X 1 vector
result_set=zeros(40);
result_nc_set=zeros(40);

%%%hK
for k=1:M
    h(:,k)=sqrt(1-sigma_e(k)^2)*(randn(Nt,1)+1i*randn(Nt,1))/sqrt(2);%Rayleigh h1
end
Pt=10^(SNR/10);


[A,B,C,D]=cal_ABCD(Nt,M,Pt,N0,h,sigma_e);


[Xc,PK,result_set]=cal_X_ZF(Nt,M,Pt,h,N0,sigma_e);



pk_zf=h*inv(h'*h);
for k=1:M
    mag(k)=norm(pk_zf(:,k));
end


pk=reshape(pk_zf.*(sqrt(PK)./mag')',[],1);


[p_max]=find_p_ZF(M,Pt,A,B,C,D,Xc,pk);

[GMI]=cal_GMI(M,A,B,C,D,p_max);


result_set(length(result_set))/log(2)
GMI







