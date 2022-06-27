% clear all


Nt=2;%nuber of transmitter antenna
K=2;%number of user
M=1;

SNR=[5, 10, 30];
%Pt=10^(SNR/10);%transmit signal power

sigma_e1=sqrt(0.1);
beta=(sigma_e1)^2;

sigma_e=randn(K,1)*0+0;%M X 1 vector
result_set_e=zeros(40);
result_nc_set_e=zeros(40);
result_zf_set_e=zeros(40);
%gam=0;
K=2;


sample=20;
for i=1:sample
    
    i
    
    
    for j=1:1:length(SNR)
               
        
    clearvars A B C D X
        
        K=2;
        Nt=K;
        
        
        H_h1=(randn(Nt,K)+1i*randn(Nt,K))/sqrt(2);%Rayleigh h1
        
        H_h=sqrt(1-beta)*H_h1;

        
        Pt=10^(SNR(j)/10);
        
        H_h=sqrt(1-beta)*H_h1;

        N0=1+beta*Pt;%gaussina noise
        N02=1;%without channel error
        
sigma_e=randn(K,1)*0+0;%M X 1 vector
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Optimization GMI
        
        [A,B,C,D]=cal_ABCD(Nt,K,Pt,N0,H_h,sigma_e*0);
        [X,result(:,j,i),converg(:,j,i)]=cal_X_converge(Nt,K,Pt,A,B,C,D);
        
        %reult/log(2)==GMI
        [GMI_X(i,j)]=cal_GMI_withX(K,A,B,C,D,X);
%         [p_max(:,i,j)]=find_p(K,Pt,A,B,C,D,X);norm(p_max(:,i,j))^2;
%         [GMI_L(i,j)]=cal_GMI(K,A,B,C,D,p_max(:,i,j));
        
    end
    
end
hold on
iter=[0,1:10];
plot(iter,mean(result(:,1,:),3),'g-*','LineWidth',1.5,'MarkerSize',8);
plot(iter,mean(result(:,2,:),3),'b-*','LineWidth',1.5,'MarkerSize',8);
plot(iter,mean(result(:,3,:),3),'r-*','LineWidth',1.5,'MarkerSize',8);
legend('SNR=5','SNR=10','SNR=20')
xlim([1 10])
ylim([min(result(:,1,1))-0.5 max(result(:,3))+0.5])
grid
xlabel('Number of iteration')
ylabel('The value of object funtion in (P3)')


