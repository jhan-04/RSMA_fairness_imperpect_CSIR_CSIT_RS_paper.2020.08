clear all
clc

Nt=2;%nuber of transmitter antenna
K=2;%number of user


SNR=[5:6:35];
%Pt=10^(SNR/10);%transmit signal power
error=sqrt(0.05);
sigma_e=error+zeros(K,1);
beta=(error)^2;
N0=1;

%gam=0;

sample=100;
for i=1:sample
    
    i
    %user1 strong channel
    %user2 week channel
    
        H_h1=(randn(Nt,K)+1i*randn(Nt,K))/sqrt(2);%Rayleigh h1
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
        
        [A,B,C,D]=cal_ABCD(Nt,K,Pt,N0,H_h,sigma_e);        
        [AA,BB,CC,DD]=cal_ABCD(Nt,K,Pt,N0,H_h,sigma_e*0);
        rho=1-abs(h1(:,i)'/norm(h1(:,i))*h2(:,i)/norm(h2(:,i)))^2;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Optimization GMI
        
        [RRR_RS0(i,j,:),RRR_sdma0(i,j,:),RRR_zf0(i,j,:),RRR_mrt0(i,j,:),RRR_noma0(i,j,:)]=CAL_SR_poroposed(Nt,K,H_h,Pt,N0,sigma_e);
        
        %%%%%%NOMA_two user
        [P_noma,RRR_noma02(i,j,:)]=GMI_NOMA_two_user(Nt,K,H_h,Pt,N0,sigma_e);
        %%%% ================================no information about channel error=
        [RRR_RS_e0(i,j,:),RRR_sdma_e0(i,j,:),RRR_zf_e0(i,j,:),RRR_mrt_e0(i,j,:),RRR_noma_e0(i,j,:)]=CAL_SR_poroposed_noinfo(Nt,K,H_h,Pt,N0,sigma_e);
        %%%%%%NOMA_two user
        [P_noma,RRR_noma]=GMI_NOMA_two_user(Nt,K,H_h,Pt,N0,sigma_e*0);
        [GMI,GMI_c,GMI_p]=cal_GMI(K,A,B,C,D,P_noma);
        RRR_noma_e02(i,j,:)=[GMI,GMI_c,GMI_p]';
        
        
                [A,B,C,D]=cal_ABCD(Nt,K,Pt,N0,H_h,sigma_e);
                p_OMA1=sqrt(Pt).*[h1(:,i);0 ;0;0 ;0]./norm(h1(:,i));
                [GMI_oma1,GMI_c_oma1,GMI_p_oma1]=cal_GMI(K,A,B,C,D, p_OMA1);
        
                p_OMA2=sqrt(Pt)*[0; 0;h2(:,i);0 ;0]./norm(h2(:,i));
                [GMI_oma2,GMI_c_oma2,GMI_p_oma2]=cal_GMI(K,A,B,C,D, p_OMA2);
                OMA_mrt(i,j,:)=[1/2.*(GMI_oma1+GMI_oma2),0,GMI_oma1,GMI_oma2];%two - user
        

    end
                    norm(h1(:,i))/norm(h2(:,i))
end
figure
hold on

plot(SNR,mean(RRR_RS0(:,:,1)),'-sr','LineWidth',2,'MarkerSize',8);plot(SNR,mean(RRR_RS_e0(:,:,1)),'--r>','LineWidth',1.5,'MarkerSize',8);
plot(SNR,mean(RRR_zf0(:,:,1)),'k:s','LineWidth',1.5,'MarkerSize',7);plot(SNR,mean(RRR_zf_e0(:,:,1)),'k-.<','LineWidth',1.5,'MarkerSize',6);
plot(SNR,mean(RRR_mrt0(:,:,1)),'g-d','LineWidth',1.5,'MarkerSize',8);plot(SNR,mean(RRR_mrt_e0(:,:,1)),'g--<','LineWidth',1.5,'MarkerSize',8);
plot(SNR,mean(OMA_mrt(:,:,1)),'c--*','LineWidth',1.5,'MarkerSize',8);
legend('RS',' RS (\sigma_e=0)','RS-ZF','RS-ZF (\sigma_e=0)',...
    'RS-MRT','RS-MRT (\sigma_e=0)','OMA-MRT')
xlabel('SNR (dB)')
ylabel('rate (bits/s/Hz)')
title({['Nt=',num2str(Nt),', M=',num2str(K),', \sigma_e=',num2str(sigma_e1),', sample=',num2str(sample)], ''})


figure
hold on
plot(SNR,mean(RRR_RS0(:,:,1)),'-sr','LineWidth',2,'MarkerSize',8);plot(SNR,mean(RRR_RS_e0(:,:,1)),'--r>','LineWidth',1.5,'MarkerSize',8);
plot(SNR,mean(RRR_sdma0(:,:,1)),'b-.s','LineWidth',1.5,'MarkerSize',8);plot(SNR,mean(RRR_sdma_e0(:,:,1)),'b--^','LineWidth',1.5,'MarkerSize',8);
plot(SNR,mean(RRR_noma0(:,:,1)),'k--s','LineWidth',1.5,'MarkerSize',8);plot(SNR,mean(RRR_noma_e0(:,:,1)),'k--<','LineWidth',1.5,'MarkerSize',8);
plot(SNR,mean(OMA_mrt(:,:,1)),'c--*','LineWidth',1.5,'MarkerSize',8);
legend('RS',' RS (\sigma_e=0)','SDMA','SDMA (\sigma_e=0)','NOMA','NOMA (\sigma_e=0)','OMA-MRT')
xlabel('SNR (dB)')
ylabel('Sum-rate (bits/s/Hz)')
title({['Nt=',num2str(Nt),', M=',num2str(K),', \sigma_e=',num2str(sigma_e1),', sample=',num2str(sample)], ''})

