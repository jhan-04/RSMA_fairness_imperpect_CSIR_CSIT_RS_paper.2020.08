%clear all


Nt=2;%nuber of transmitter antenna
K=2;%number of user
M=1;

SNR=20;
%Pt=10^(SNR/10);%transmit signal power

sigma_e1=0.1;
beta=(sigma_e1)^2;

sigma_e=randn(K,1)*0+0;%M X 1 vector
result_set_e=zeros(40);
result_nc_set_e=zeros(40);
result_zf_set_e=zeros(40);
%gam=0;
KK=[2,4,8];


sample=1;
for i=1:sample
    
    i
    
    
    for j=1:1:length(KK)
               
        
    clearvars A B C D X
        
        K=KK(j)
        Nt=K;
        
        
        H_h1=(randn(Nt,K)+1i*randn(Nt,K))/sqrt(2);%Rayleigh h1
        
        H_h=sqrt(1-beta)*H_h1;

        
        Pt=10^(SNR/10);
        
        H_h=sqrt(1-beta)*H_h1;

        N0=1+beta*Pt;%gaussina noise
        N02=1;%without channel error
        
sigma_e=randn(K,1)*0+0;%M X 1 vector
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Optimization GMI
        
        [A,B,C,D]=cal_ABCD(Nt,K,Pt,N0,H_h,sigma_e*0);
        [X,result(:,j)]=cal_X_converge(Nt,K,Pt,A,B,C,D);
        
        %reult/log(2)==GMI
        [GMI_X(i,j)]=cal_GMI_withX(K,A,B,C,D,X);
%         [p_max(:,i,j)]=find_p(K,Pt,A,B,C,D,X);norm(p_max(:,i,j))^2;
%         [GMI_L(i,j)]=cal_GMI(K,A,B,C,D,p_max(:,i,j));
        
    end
    
end
hold on
iter=[0,1:20];
plot(iter,result(:,1),'g-*','LineWidth',1.5,'MarkerSize',8);
plot(iter,result(:,2),'b-*','LineWidth',1.5,'MarkerSize',8);
plot(iter,result(:,3),'r-*','LineWidth',1.5,'MarkerSize',8);
legend('K,N_t=2','K,N_t=4','K,N_t=8')
xlim([1 20])
ylim([2 max(result(:,3))+2])
grid