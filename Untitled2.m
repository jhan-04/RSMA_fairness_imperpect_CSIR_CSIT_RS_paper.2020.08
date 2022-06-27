% clear all
% clc

Nt=3;%nuber of transmitter antenna
K=3;%number of user


SNR=[5:6:35];
%Pt=10^(SNR/10);%transmit signal power
error=(0.15);
sigma_e=error+zeros(K,1);
beta=(error)^2;
N0=1;

%gam=0;

sample=20;
for i=11:sample
    
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
        [A_t,B_t,C_t,D_t]=cal_ABCD(Nt,K,Pt,N0,H_h,zeros(K,1));
        [X,CC,result]=cal_X(Nt,K,Pt,A,B,C,D);
        %result_set(1:length(result))=result;
        [p_max(:,i,j)]=find_p(K,Pt,A,B,C,D,X,CC);%norm(p_max(:,i,j))^2;
        
        [GMI_min,GMI_c,GMI_p]=cal_GMI_min(K,A,B,C,D,p_max(:,i,j),CC);
        
        RRR_RS(i,j,:)=[GMI_min,GMI_c,GMI_p]';
        CCC(i,j,:)=CC;
        
        
        
        %%%%%imperfect CSIT
         [GMI_r,GMI_c_r,GMI_p_r]=cal_GMI_min(K,A_t,B_t,C_t,D_t,p_max(:,i,j),CC);
        RRR_RS_imperfect_CSIT(i,j,:)=[GMI_r,GMI_c_r,GMI_p_r]';
        
        
        
        %%SDMA
        [A,B,C,D]=cal_ABCD(Nt,K,Pt,N0,H_h,sigma_e);
        [X_nc,result_nc]=cal_X_no_common(Nt,K,Pt,A,B,C,D);
        [p_nc(:,i,j)]=find_p(K,Pt,A,B,C,D,X_nc,CC);%norm(p_nc)^2;
        [GMI_min_SDMA,GMI_c_SDMA,GMI_p_SDMA]=cal_GMI_min(K,A,B,C,D,p_nc(:,i,j),rand(K,1));
        
        
        RRR_sdma(i,j,:)= [GMI_min_SDMA,GMI_c_SDMA,GMI_p_SDMA]';
        
        [GMI_r_sdma,GMI_c_r_sdma,GMI_p_r_sdma]=cal_GMI_min(K,A_t,B_t,C_t,D_t,p_nc(:,i,j),CC);
        RRR_sdma_imperfect_CSIT(i,j,:)=[GMI_r_sdma,GMI_c_r_sdma,GMI_p_r_sdma]';
        
        %NOMA
        [A,B,C,D]=cal_ABCD(Nt,K,Pt,N0,H_h,sigma_e);
        [X_N,result_set_N]=cal_X_NOMA(Nt,K,H_h,A,B,C,D);
        [p_N(:,i,j)]=find_p(K,Pt,A,B,C,D,X_N,[0;100]);
        [GMI_min_noma,GMI_c_noma,GMI_p_noma]=cal_GMI_min(K,A,B,C,D,p_N(:,i,j),[0;100]);
        
        RRR_noma(i,j,:)=[GMI_min_noma,GMI_c_noma,GMI_p_noma]';
       
        [GMI_r_noma,GMI_c_r_noma,GMI_p_r_noma]=cal_GMI(K,A_t,B_t,C_t,D_t,p_N(:,i,j));
        RRR_noma_imperfect_CSIT(i,j,:)=[GMI_r_noma,GMI_c_r_noma,GMI_p_r_noma]';
        
    end
    %norm(h1(:,i))/norm(h2(:,i))
end
figure
hold on

plot(SNR,mean(RRR_RS(:,:,1)),'-sr','LineWidth',2,'MarkerSize',8);%plot(SNR,mean(RRR_RS_e0(:,:,1)),'--r>','LineWidth',1.5,'MarkerSize',8);
plot(SNR,mean(RRR_sdma(:,:,1)),'-sb','LineWidth',2,'MarkerSize',8);
plot(SNR,mean(RRR_noma(:,:,1)),'-sg','LineWidth',2,'MarkerSize',8);

plot(SNR,mean(RRR_RS(:,:,3)),'--xr','LineWidth',2,'MarkerSize',6);%plot(SNR,mean(RRR_RS_e0(:,:,1)),'--r>','LineWidth',1.5,'MarkerSize',8);
plot(SNR,mean(RRR_sdma(:,:,3)),'--xb','LineWidth',2,'MarkerSize',6);
plot(SNR,mean(RRR_noma(:,:,3)),'--xg','LineWidth',2,'MarkerSize',6);

plot(SNR,mean(RRR_RS(:,:,4)),'--or','LineWidth',2,'MarkerSize',5);%plot(SNR,mean(RRR_RS_e0(:,:,1)),'--r>','LineWidth',1.5,'MarkerSize',8);
plot(SNR,mean(RRR_sdma(:,:,4)),'--ob','LineWidth',2,'MarkerSize',5);
plot(SNR,mean(RRR_noma(:,:,4)),'--og','LineWidth',2,'MarkerSize',5);


plot(SNR,mean(RRR_RS(:,:,2).*CCC(:,:,1)./sum(CCC(:,:,:),3)),':^r','LineWidth',2,'MarkerSize',4);%plot(SNR,mean(RRR_RS_e0(:,:,1)),'--r>','LineWidth',1.5,'MarkerSize',8);
plot(SNR,mean(RRR_sdma(:,:,2).*CCC(:,:,1)./sum(CCC(:,:,:),3)),':^b','LineWidth',2,'MarkerSize',4);
plot(SNR,mean(RRR_noma(:,:,2).*CCC(:,:,1)./sum(CCC(:,:,:),3)),':^g','LineWidth',2,'MarkerSize',4);

plot(SNR,mean(RRR_RS(:,:,2).*CCC(:,:,2)./sum(CCC(:,:,:),3)),':>r','LineWidth',2,'MarkerSize',4);%plot(SNR,mean(RRR_RS_e0(:,:,1)),'--r>','LineWidth',1.5,'MarkerSize',8);
plot(SNR,mean(RRR_sdma(:,:,2).*CCC(:,:,2)./sum(CCC(:,:,:),3)),':>b','LineWidth',2,'MarkerSize',4);
plot(SNR,mean(RRR_noma(:,:,2).*CCC(:,:,2)./sum(CCC(:,:,:),3)),':>g','LineWidth',2,'MarkerSize',4);



legend('RS','SDMA','NOMA','RS-pivate-1','SDMA-pivate-1','NOMA-pivate-1','RS-pivate-2','SDMA-pivate-2','NOMA-pivate-2','RS-common-1','SDMA-common-1','NOMA-common-1','RS-common-2','SDMA-common-2','NOMA-common-2')
xlabel('SNR (dB)')
ylabel('min rate (bits/s/Hz)')
title({['Nt=',num2str(Nt),', K=',num2str(K),', \sigma_e ^2=',num2str(error^2),', sample=',num2str(sample)], 'imperfect CSIR'})




figure
hold on

plot(SNR,mean(RRR_RS_imperfect_CSIT(:,:,1)),'-sr','LineWidth',2,'MarkerSize',8);%plot(SNR,mean(RRR_RS_e0(:,:,1)),'--r>','LineWidth',1.5,'MarkerSize',8);
plot(SNR,mean(RRR_sdma_imperfect_CSIT(:,:,1)),'-sb','LineWidth',2,'MarkerSize',8);
%plot(SNR,mean(RRR_noma_imperfect_CSIT(:,:,1)),'-sg','LineWidth',2,'MarkerSize',8);

legend('RS','SDMA')%,'NOMA')
xlabel('SNR (dB)')
ylabel('min rate (bits/s/Hz)')
title({['Nt=',num2str(Nt),', K=',num2str(K),', \sigma_e ^2=',num2str(error^2),', sample=',num2str(sample)], 'imperfect CSIT'})




