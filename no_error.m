clear all
clc

Nt=2;%nuber of transmitter antenna
K=2;%number of user


SNR=[5:6:35];
%Pt=10^(SNR/10);%transmit signal power
error=sqrt(0);
sigma_e=error+zeros(K,1);
beta=(error)^2;
N0=1;

%gam=0;

sample=50;
for i=1:sample
    
    i
    %user1 strong channel
    %user2 week channel
    
    H_h1=(randn(Nt,K)+1i*randn(Nt,K))/sqrt(2);%.*[1;2];%Rayleigh h1
    H_h=sqrt(1-beta)*H_h1;
    if norm(H_h(:,1))>=norm(H_h(:,2))
        h1(:,i)=H_h(:,1); h2(:,i)=H_h(:,2);
    else
        h1(:,i)=H_h(:,2); h2(:,i)=H_h(:,1);
    end
    H_h(:,1)=h1(:,i);
    H_h(:,2)=h2(:,i);
    
    for j=1:1:length(SNR)
        
        
        snr=SNR(j)
        Pt=10^(snr/10);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Optimization GMI
        
        

        %%RSMA
        [A,B,C,D]=cal_ABCD(Nt,K,Pt,N0,H_h,sigma_e);
        [X,CC,result]=cal_X(Nt,K,Pt,A,B,C,D);
        %result_set(1:length(result))=result;
        [p_max(:,i,j)]=find_p(K,Pt,A,B,C,D,X,CC);%norm(p_max(:,i,j))^2;
        
        [GMI_min,GMI_c,GMI_p]=cal_GMI_min(K,A,B,C,D,p_max(:,i,j),CC);
        
        RRR_RS(i,j,:)=[GMI_min,GMI_c,GMI_p]';
        
        CCC(i,j,:)=CC;
        
        
        
        
        
        
        
        
        %%SDMA
        [A,B,C,D]=cal_ABCD(Nt,K,Pt,N0,H_h,sigma_e);
        [X_nc,result_nc]=cal_X_no_common(Nt,K,Pt,A,B,C,D);
        [p_nc(:,i,j)]=find_p(K,Pt,A,B,C,D,X_nc,CC);%norm(p_nc)^2;
        [GMI_min_SDMA,GMI_c_SDMA,GMI_p_SDMA]=cal_GMI_min(K,A,B,C,D,p_nc(:,i,j),rand(K,1));
        
        
        RRR_sdma(i,j,:)= [GMI_min_SDMA,GMI_c_SDMA,GMI_p_SDMA]';
        
        %NOMA
        [A,B,C,D]=cal_ABCD(Nt,K,Pt,N0,H_h,sigma_e);
        [X_N,result_set_N]=cal_X_NOMA(Nt,K,H_h,A,B,C,D);
        [p_N(:,i,j)]=find_p(K,Pt,A,B,C,D,X_N,[0;100]);
        [GMI_min_noma,GMI_c_noma,GMI_p_noma]=cal_GMI_min(K,A,B,C,D,p_N(:,i,j),[0;100]);
        
        RRR_noma(i,j,:)=[GMI_min_noma,GMI_c_noma,GMI_p_noma]';
        
        
    end
    %norm(h1(:,i))/norm(h2(:,i))
end
figure
hold on

plot(SNR,mean(RRR_RS(:,:,1)),'-sr','LineWidth',2,'MarkerSize',8);%plot(SNR,mean(RRR_RS_e0(:,:,1)),'--r>','LineWidth',1.5,'MarkerSize',8);
plot(SNR,mean(RRR_sdma(:,:,1)),'-sb','LineWidth',2,'MarkerSize',8);
plot(SNR,mean(RRR_noma(:,:,1)),'-sg','LineWidth',2,'MarkerSize',8);