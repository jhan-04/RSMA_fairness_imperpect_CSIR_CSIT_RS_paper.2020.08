% clear all


Nt=2;%nuber of transmitter antenna
K=2;%number of user
M=1;

SNR=40:10:60;
%Pt=10^(SNR/10);%transmit signal power

sigma_e1=0.1;
beta=(sigma_e1)^2;

sigma_e=randn(K,1)*0+0;%M X 1 vector
result_set_e=zeros(40);
result_nc_set_e=zeros(40);
result_zf_set_e=zeros(40);
%gam=0;

sample=3;
for i=1:sample
    
    i
    
    H_h1=(randn(Nt,K)+1i*randn(Nt,K))/sqrt(2);%Rayleigh h1
    
    H_h=sqrt(1-beta)*H_h1;
    if norm(H_h(:,1))>=norm(H_h(:,2))
        h1(:,i)=H_h(:,1); h2(:,i)=H_h(:,2);
    else
        h1(:,i)=H_h(:,2); h2(:,i)=H_h(:,1);
    end
    
    
    %     Rho(i)=rho;
    
    for j=1:1:length(SNR)
        
        
        snr=SNR(j)
        Pt=10^(snr/10);
        
        H_h=sqrt(1-beta)*H_h1;
        
        
        
        rho=1-abs(h1(:,i)'/norm(h1(:,i))*h2(:,i)/norm(h2(:,i)))^2;
        N0=1+beta*Pt;%gaussina noise
        N02=1;%without channel error
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Optimization GMI
        [RRR_G(i,j,:),RRR_sdma(i,j,:),RRR_zf(i,j,:),RRR_mrt(i,j,:)]=CAL_SR_poroposed_multiuser(Nt,K,H_h,H_h,Pt,M,N0,sigma_e*0);
        %================================no information about channel error=
        [RRR_G_e(i,j,:),RRR_sdma_e(i,j,:),RRR_zf_e(i,j,:),RRR_mrt_e(i,j,:)]=CAL_SR_poroposed_noinfo_multiuser(Nt,K,H_h,H_h,Pt,M,N0,N02,sigma_e*0);


        
        %%MMSE
        [pc_1,pk_1]=Alg_1(Nt,K,H_h,H_h,Pt,M,N0); ESR(i,j)=Alg_1_cal_ESR(pc_1,pk_1,H_h,0,M,K,N0);%paper
        [pc_con,pk_con,Bound(i,j)]=Alg_1_con(Nt,K,H_h,Pt,N0,sigma_e(1)^2*0); ESR_con(i,j)=Alg_1_cal_ESR(pc_con,pk_con,H_h,0,M,K,N0);%paper
        
        [pc_noRS,pk_noRS]=Alg_1_noRS2(Nt,K,H_h,H_h,Pt,M,N0);
        ESR_noRS(i,j)=Alg_1_cal_ESR(pc_noRS,pk_noRS,H_h,0,M,K,N0);%paper
        
    end
    
end
figure(1)
hold on
%ahevable rate
%RS imperfect with error info&%imperfect with no error info
plot(SNR,mean(RRR_G(:,:,1)),'-sr','LineWidth',2,'MarkerSize',8);plot(SNR,mean(RRR_G_e(:,:,1)),'--r>','LineWidth',1.5,'MarkerSize',8);

plot(SNR,mean(RRR_sdma(:,:,1)),'b-.s','LineWidth',1.5,'MarkerSize',8);plot(SNR,mean(RRR_sdma_e(:,:,1)),'b--^','LineWidth',1.5,'MarkerSize',8);

plot(SNR,mean(RRR_zf(:,:,1)),'k:s','LineWidth',1.5,'MarkerSize',7);plot(SNR,mean(RRR_zf_e(:,:,1)),'k-.<','LineWidth',1.5,'MarkerSize',6);

plot(SNR,mean(RRR_mrt(:,:,1)),'g-d','LineWidth',1.5,'MarkerSize',8);plot(SNR,mean(RRR_mrt_e(:,:,1)),'g--<','LineWidth',1.5,'MarkerSize',8);

%plot(SNR,mean(RRR_noma(:,:,1)),'m--d','LineWidth',1.5,'MarkerSize',8);plot(SNR,mean(RRR_noma_e(:,:,1)),'m--<','LineWidth',1.5,'MarkerSize',8);




legend('RS',' RS_{no error info}','SDMA','SDMA_{no error info}','RS-ZF','RS-ZF_{no error info}',...
    'RS-MRT','RS-MRT_{no error info}')

xlabel('SNR (dB)')
ylabel('rate (bits/s/Hz)')
title(['Nt=',num2str(Nt),', M=',num2str(K),', \sigma_e=',num2str(sigma_e1),', sample=',num2str(sample)])

hold on
%MMSE
%plot(SNR,mean(ESR),'m')
plot(SNR,mean(ESR_con),'m--o')
