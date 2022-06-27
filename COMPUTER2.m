%clear all


Nt=2;%nuber of transmitter antenna
K=2;%number of user
M=1;

SNR=5:6:35;
%Pt=10^(SNR/10);%transmit signal power

sigma_e1=sqrt(0.1);
beta=(sigma_e1)^2;

sigma_e=randn(K,1)*0+0;%M X 1 vector
result_set_e=zeros(40);
result_nc_set_e=zeros(40);
result_zf_set_e=zeros(40);
%gam=0;

sample=19;
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
        
        
        
        
        rho=1-abs(h1(:,i)'/norm(h1(:,i))*h2(:,i)/norm(h2(:,i)))^2;
        N0=1+beta*Pt;%gaussina noise
        N02=1;%without channel error
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Optimization GMI
        %[RRR_G(i,j,:),RRR_sdma(i,j,:),RRR_zf(i,j,:),RRR_mrt(i,j,:),RRR_noma(i,j,:)]=CAL_SR_poroposed(Nt,K,H_h,H_h,Pt,M,N0,sigma_e*0);
        %================================no information about channel error=
        [RRR_G_e(i,j,:),RRR_sdma_e(i,j,:),RRR_zf_e(i,j,:),RRR_mrt_e(i,j,:),RRR_noma_e(i,j,:)]=CAL_SR_poroposed_noinfo(Nt,K,H_h,H_h,Pt,M,N0,N02,sigma_e*0);
        %%%%Zero forcing
%         [RRR_zf_i(i,j,:),RRR_zf_p(i,j,:),RRR_zf_ni(i,j,:),RRR_zf_np(i,j,:),RRR_zf_si(i,j,:),RRR_zf_sp(i,j,:),...
%             RRR_zf_mi(i,j,:),RRR_zf_mp(i,j,:),RRR_mrt1(i,j,:),RRR_mrt2(i,j,:)]=CAL_SR_ZF_poroposed_twoUser(H_h,Pt,beta);
%         
% RRR_noma(i,j,:);
%  RRR_mrt1(i,j,:);       
%         %%MMSE
%         [pc_1,pk_1]=Alg_1(Nt,K,H_h,H_h,Pt,M,N0); ESR(i,j)=Alg_1_cal_ESR(pc_1,pk_1,H_h,0,M,K,N0);%paper
%         [pc_con,pk_con,Bound(i,j)]=Alg_1_con(Nt,K,H_h,Pt,N0,sigma_e(1)^2*0); ESR_con(i,j)=Alg_1_cal_ESR(pc_con,pk_con,H_h,0,M,K,N0);%paper
%         
%         [pc_noRS,pk_noRS]=Alg_1_noRS2(Nt,K,H_h,H_h,Pt,M,N0);
%         ESR_noRS(i,j)=Alg_1_cal_ESR(pc_noRS,pk_noRS,H_h,0,M,K,N0);%paper
%         
    end
    
end
figure(1)
hold on
%ahevable rate
%RS imperfect with error info&%imperfect with no error info
plot(SNR,mean(RRR_G(:,:,1)),'-sr','LineWidth',2,'MarkerSize',8);plot(SNR,mean(RRR_G_e(:,:,1)),'--r>','LineWidth',1.5,'MarkerSize',8);

plot(SNR,mean(RRR_sdma(:,:,1)),'b-.s','LineWidth',1.5,'MarkerSize',8);plot(SNR,mean(RRR_sdma_e(:,:,1)),'b--^','LineWidth',1.5,'MarkerSize',8);

plot(SNR,mean(RRR_zf(:,:,1)),'k:s','LineWidth',1.5,'MarkerSize',7);plot(SNR,mean(RRR_zf_e(:,:,1)),'k-.<','LineWidth',1.5,'MarkerSize',6);

plot(SNR,mean(RRR_mrt(:,:,1)),'y-d','LineWidth',1.5,'MarkerSize',8);plot(SNR,mean(RRR_mrt_e(:,:,1)),'y--<','LineWidth',1.5,'MarkerSize',8);

plot(SNR,mean(RRR_noma(:,:,1)),'m--d','LineWidth',1.5,'MarkerSize',8);plot(SNR,mean(RRR_noma_e(:,:,1)),'m--<','LineWidth',1.5,'MarkerSize',8);



%%%%ZeroForcing
%%RS
plot(SNR,mean(RRR_zf_i(:,:,1)),'k','LineWidth',1.5,'MarkerSize',4);plot(SNR,mean(RRR_zf_p(:,:,1)),'k--o','LineWidth',1.5,'MarkerSize',8);
%%NOMA
plot(SNR,mean(RRR_zf_ni(:,:,1)),'g','LineWidth',1.5,'MarkerSize',8);plot(SNR,mean(RRR_zf_np(:,:,1)),'g--o','LineWidth',1.5,'MarkerSize',8);
%SDMA
plot(SNR,mean(RRR_zf_si(:,:,1)),'b','LineWidth',1.5,'MarkerSize',8);plot(SNR,mean(RRR_zf_sp(:,:,1)),'b--o','LineWidth',1.5,'MarkerSize',8);
%MRT
plot(SNR,mean(RRR_mrt1(:,:,1)),'c--*','LineWidth',1.5,'MarkerSize',8);plot(SNR,mean(RRR_mrt2(:,:,1)),'c--+','LineWidth',1.5,'MarkerSize',8);
%multicasting
plot(SNR,mean(RRR_zf_mi(:,:,1)),'m');plot(SNR,mean(RRR_zf_mp(:,:,1)),'m--o');

legend('RS',' RS_{no error info}','SDMA','SDMA_{no error info}','RS-ZF','RS-ZF_{no error info}',...
    'RS-MRT','RS-MRT_{no error info}','NOMA','NOMA_{no error info}',...
    'RS-ZF','RS-ZF_{no error info}','NOMA-ZF','NOMA-ZF_{no error info}','SDMA-ZF','SDMA-ZF_{no error info}',...
    'OMA-MRT1','OMA-MRT2','multicasting','multicasting_{no error info}')

xlabel('SNR (dB)')
ylabel('rate (bits/s/Hz)')
title(['Nt=',num2str(Nt),', M=',num2str(K),', \sigma_e=',num2str(sigma_e1),', sample=',num2str(sample)])

hold on
%MMSE
plot(SNR,mean(ESR),'m')
plot(SNR,mean(ESR_con),'m--o')

figure(2)
hold on
plot(SNR,mean(RRR_G(:,:,1)),'-sr','LineWidth',2,'MarkerSize',8);plot(SNR,mean(RRR_G_e(:,:,1)),'--r>','LineWidth',1.5,'MarkerSize',8);
plot(SNR,mean(RRR_sdma(:,:,1)),'b-.s','LineWidth',1.5,'MarkerSize',8);plot(SNR,mean(RRR_sdma_e(:,:,1)),'b--^','LineWidth',1.5,'MarkerSize',8);
plot(SNR,mean(RRR_noma(:,:,1)),'g--d','LineWidth',1.5,'MarkerSize',8);plot(SNR,mean(RRR_noma_e(:,:,1)),'g--<','LineWidth',1.5,'MarkerSize',8);
plot(SNR,mean(RRR_zf_mi(:,:,1)),'m-o','LineWidth',1.5,'MarkerSize',8);
plot(SNR,mean(RRR_mrt1(:,:,1)),'c--*','LineWidth',1.5,'MarkerSize',8);
legend('RS',' RS_{no error info}','SDMA','SDMA_{no error info}','NOMA','NOMA_{no error info}',...
    'multicasting','OMA-MRT')
xlabel('SNR (dB)')
ylabel('rate (bits/s/Hz)')
title({['Nt=',num2str(Nt),', M=',num2str(K),', \sigma_e=',num2str(sigma_e1),', sample=',num2str(sample)], ['Compared with other multiple access strategies']})



figure(3)
hold on
plot(SNR,mean(RRR_G(:,:,1)),'-sr','LineWidth',2,'MarkerSize',8);plot(SNR,mean(RRR_G_e(:,:,1)),'--r>','LineWidth',1.5,'MarkerSize',8);
plot(SNR,mean(RRR_zf(:,:,1)),'k:s','LineWidth',1.5,'MarkerSize',7);plot(SNR,mean(RRR_zf_e(:,:,1)),'k-.<','LineWidth',1.5,'MarkerSize',6);
plot(SNR,mean(RRR_mrt(:,:,1)),'g-d','LineWidth',1.5,'MarkerSize',8);plot(SNR,mean(RRR_mrt_e(:,:,1)),'g--<','LineWidth',1.5,'MarkerSize',8);
plot(SNR,mean(RRR_zf_si(:,:,1)),'b','LineWidth',1.5,'MarkerSize',8);plot(SNR,mean(RRR_zf_sp(:,:,1)),'b--o','LineWidth',1.5,'MarkerSize',8);
plot(SNR,mean(RRR_mrt1(:,:,1)),'c--*','LineWidth',1.5,'MarkerSize',8);
legend('RS',' RS_{no error info}','RS-ZF','RS-ZF_{no error info}',...
    'RS-MRT','RS-MRT_{no error info}','SDMA-ZF','SDMA-ZF_{no error info}','OMA-MRT')
xlabel('SNR (dB)')
ylabel('rate (bits/s/Hz)')
title({['Nt=',num2str(Nt),', M=',num2str(K),', \sigma_e=',num2str(sigma_e1),', sample=',num2str(sample)], 'Compared with other precoding vector'})

