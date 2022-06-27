% clear all


Nt=2;%nuber of transmitter antenna
K=2;%number of user
M=1;
SNR=30;
Pt=10^(SNR/10);
ERROR=([0.01,0.02,0.05,0.1,0.2]);
%Pt=10^(SNR/10);%transmit signal power


sigma_e=randn(K,1)*0+0;%M X 1 vector
result_set_e=zeros(40);
result_nc_set_e=zeros(40);
result_zf_set_e=zeros(40);
%gam=0;

sample=100;
for i=11:sample
    
    i
    
    H_h1=(randn(Nt,K)+1i*randn(Nt,K))/sqrt(2);%Rayleigh h1
    
    
    
    %     Rho(i)=rho;
    
    for j=1:1:length(ERROR)
        sigma_e1=sqrt(ERROR(j))
        beta=(sigma_e1)^2;
        
        H_h=sqrt(1-beta)*H_h1;
        
        if norm(H_h(:,1))>=norm(H_h(:,2))
            h1(:,i)=H_h(:,1); h2(:,i)=H_h(:,2);
        else
            h1(:,i)=H_h(:,2); h2(:,i)=H_h(:,1);
        end
               H_h(:,1)=h1(:,i);
        H_h(:,2)=h2(:,i);
        
        
        
        rho=1-abs(h1(:,i)'/norm(h1(:,i))*h2(:,i)/norm(h2(:,i)))^2;
        N0=1+beta*Pt;%gaussina noise
        N02=1;%without channel error
               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Optimization GMI
        [RRR_G(i,j,:),RRR_sdma(i,j,:),RRR_zf(i,j,:),RRR_mrt(i,j,:),RRR_noma(i,j,:)]=CAL_SR_poroposed(Nt,K,H_h,H_h,Pt,M,N0,sigma_e*0);
        %%%% ================================no information about channel error=
        [RRR_G_e(i,j,:),RRR_sdma_e(i,j,:),RRR_zf_e(i,j,:),RRR_mrt_e(i,j,:),RRR_noma_e(i,j,:)]=CAL_SR_poroposed_noinfo(Nt,K,H_h,H_h,Pt,M,N0,N02,sigma_e*0);
        %%Zero forcing
        [RRR_zf_i(i,j,:),RRR_zf_p(i,j,:),RRR_zf_ni(i,j,:),RRR_zf_np(i,j,:),RRR_zf_si(i,j,:),RRR_zf_sp(i,j,:),...
            RRR_zf_mi(i,j,:),RRR_zf_mp(i,j,:),RRR_mrt1(i,j,:),RRR_mrt2(i,j,:)]=CAL_SR_ZF_poroposed_twoUser(H_h,Pt,beta);
       

        %%%%%NOMA_two user
        
        [A,B,C,D]=cal_ABCD(Nt,K,Pt,N0,H_h,sigma_e*0);
        [X_N,result_set_N]=cal_X_NOMA2(Nt,K,H_h,A,B,C,D);
        [p_N(:,i,j)]=find_p_NOMA(K,Pt,A,B,C,D,X_N,H_h);
        [RRR_noma2(i,j,:)]=cal_GMI(K,A,B,C,D,p_N(:,i,j));

        
        
        [AA,BB,CC,DD]=cal_ABCD(Nt,K,Pt,N02,H_h,sigma_e*0);
        [X_N_e,result_set_N]=cal_X_NOMA2(Nt,K,H_h,AA,BB,CC,DD);
        [p_N_e(:,i,j)]=find_p_NOMA(K,Pt,AA,BB,CC,DD,X_N_e,H_h);
        [RRR_noma_e2(i,j,:)]=cal_GMI(K,A,B,C,D,p_N_e(:,i,j));
         
        
        %         %%MMSE
        %         [pc_1,pk_1]=Alg_1(Nt,K,H_h,H_h,Pt,M,N0); ESR(i,j)=Alg_1_cal_ESR(pc_1,pk_1,H_h,0,M,K,N0);%paper
        %         [pc_con,pk_con,Bound(i,j)]=Alg_1_con(Nt,K,H_h,Pt,N0,sigma_e(1)^2*0); ESR_con(i,j)=Alg_1_cal_ESR(pc_con,pk_con,H_h,0,M,K,N0);%paper
        %
        %         [pc_noRS,pk_noRS]=Alg_1_noRS2(Nt,K,H_h,H_h,Pt,M,N0);
        %         ESR_noRS(i,j)=Alg_1_cal_ESR(pc_noRS,pk_noRS,H_h,0,M,K,N0);%paper
        
    end
    
end
figure(1)
hold on
%ahevable rate
%RS imperfect with error info&%imperfect with no error info

plot(ERROR,mean(RRR_G(:,:,1)),'-sr','LineWidth',2,'MarkerSize',8);plot(ERROR,mean(RRR_G_e(:,:,1)),'--r>','LineWidth',1.5,'MarkerSize',8);

plot(ERROR,mean(RRR_sdma(:,:,1)),'b-.s','LineWidth',1.5,'MarkerSize',8);plot(ERROR,mean(RRR_sdma_e(:,:,1)),'b--^','LineWidth',1.5,'MarkerSize',8);

plot(ERROR,mean(RRR_zf(:,:,1)),'k:s','LineWidth',1.5,'MarkerSize',7);plot(ERROR,mean(RRR_zf_e(:,:,1)),'k-.<','LineWidth',1.5,'MarkerSize',6);

plot(ERROR,mean(RRR_mrt(:,:,1)),'y-d','LineWidth',1.5,'MarkerSize',8);plot(ERROR,mean(RRR_mrt_e(:,:,1)),'y--<','LineWidth',1.5,'MarkerSize',8);

plot(ERROR,mean(RRR_noma(:,:,1)),'m--d','LineWidth',1.5,'MarkerSize',8);plot(ERROR,mean(RRR_noma_e(:,:,1)),'m--<','LineWidth',1.5,'MarkerSize',8);



%%%%ZeroForcing
%%RS
plot(ERROR,mean(RRR_zf_i(:,:,1)),'k','LineWidth',1.5,'MarkerSize',4);plot(ERROR,mean(RRR_zf_p(:,:,1)),'k--o','LineWidth',1.5,'MarkerSize',8);
%%NOMA
plot(ERROR,mean(RRR_zf_ni(:,:,1)),'g','LineWidth',1.5,'MarkerSize',8);plot(ERROR,mean(RRR_zf_np(:,:,1)),'g--o','LineWidth',1.5,'MarkerSize',8);
%SDMA
plot(ERROR,mean(RRR_zf_si(:,:,1)),'b','LineWidth',1.5,'MarkerSize',8);plot(ERROR,mean(RRR_zf_sp(:,:,1)),'b--o','LineWidth',1.5,'MarkerSize',8);
%MRT
plot(ERROR,mean(RRR_mrt1(:,:,1)),'c--*','LineWidth',1.5,'MarkerSize',8);plot(ERROR,mean(RRR_mrt2(:,:,1)),'c--+','LineWidth',1.5,'MarkerSize',8);
%multicasting
plot(ERROR,mean(RRR_zf_mi(:,:,1)),'m');plot(ERROR,mean(RRR_zf_mp(:,:,1)),'m--o');

legend('RS',' RS_{no error info}','SDMA','SDMA_{no error info}','RS-ZF','RS-ZF_{no error info}',...
    'RS-MRT','RS-MRT_{no error info}','NOMA','NOMA_{no error info}',...
    'RS-ZF','RS-ZF_{no error info}','NOMA-ZF','NOMA-ZF_{no error info}','SDMA-ZF','SDMA-ZF_{no error info}',...
    'OMA-MRT1','OMA-MRT2','multicasting','multicasting_{no error info}')

xlabel('SNR (dB)')
ylabel('rate (bits/s/Hz)')
title(['Nt=',num2str(Nt),', M=',num2str(K),', \sigma_e=',num2str(sigma_e1),', sample=',num2str(sample)])

hold on
% %MMSE
% %plot(SNR,mean(ESR),'m')
% plot(ERROR,mean(ESR_con),'m--o')




figure(2)
hold on
plot(ERROR,mean(RRR_G(:,:,1)),'-sr','LineWidth',2,'MarkerSize',8);plot(ERROR,mean(RRR_G_e(:,:,1)),'--r>','LineWidth',1.5,'MarkerSize',8);
plot(ERROR,mean(RRR_sdma(:,:,1)),'b-.s','LineWidth',1.5,'MarkerSize',8);plot(ERROR,mean(RRR_sdma_e(:,:,1)),'b--^','LineWidth',1.5,'MarkerSize',8);
plot(ERROR,mean(RRR_noma(:,:,1)),'m--d','LineWidth',1.5,'MarkerSize',8);plot(ERROR,mean(RRR_noma_e(:,:,1)),'m--<','LineWidth',1.5,'MarkerSize',8);
plot(ERROR,mean(RRR_zf_mi(:,:,1)),'m');plot(ERROR,mean(RRR_zf_mp(:,:,1)),'m--o');
plot(ERROR,mean(RRR_mrt1(:,:,1)),'c--*','LineWidth',1.5,'MarkerSize',8);
legend('RS',' RS_{no error info}','SDMA','SDMA_{no error info}','NOMA','NOMA_{no error info}',...
    'multicasting','multicasting_{no error info}','OMA-MRT1')
xlabel('SNR (dB)')
ylabel('rate (bits/s/Hz)')
title({['Nt=',num2str(Nt),', M=',num2str(K),', \sigma_e=',num2str(sigma_e1),', sample=',num2str(sample)], ['Compared with other multiple access strategies']})



figure(3)
hold on
plot(ERROR,mean(RRR_G(:,:,1)),'-sr','LineWidth',2,'MarkerSize',8);plot(ERROR,mean(RRR_G_e(:,:,1)),'--r>','LineWidth',1.5,'MarkerSize',8);
plot(ERROR,mean(RRR_zf(:,:,1)),'k:s','LineWidth',1.5,'MarkerSize',7);plot(ERROR,mean(RRR_zf_e(:,:,1)),'k-.<','LineWidth',1.5,'MarkerSize',6);
plot(ERROR,mean(RRR_mrt(:,:,1)),'y-d','LineWidth',1.5,'MarkerSize',8);plot(ERROR,mean(RRR_mrt_e(:,:,1)),'y--<','LineWidth',1.5,'MarkerSize',8);
plot(ERROR,mean(RRR_zf_si(:,:,1)),'b','LineWidth',1.5,'MarkerSize',8);plot(ERROR,mean(RRR_zf_sp(:,:,1)),'b--o','LineWidth',1.5,'MarkerSize',8);
plot(ERROR,mean(RRR_mrt1(:,:,1)),'c--*','LineWidth',1.5,'MarkerSize',8);
legend('RS',' RS_{no error info}','RS-ZF','RS-ZF_{no error info}',...
    'RS-MRT','RS-MRT_{no error info}','SDMA-ZF','SDMA-ZF_{no error info}','OMA-MRT1')
xlabel('SNR (dB)')
ylabel('rate (bits/s/Hz)')
title({['Nt=',num2str(Nt),', M=',num2str(K),', \sigma_e=',num2str(sigma_e1),', sample=',num2str(sample)], 'Compared with other precoding vector'})






figure(4)
hold on
%X = categorical({num2str(ERROR(1)),num2str(ERROR(2)),num2str(ERROR(3)),num2str(ERROR(4)),num2str(ERROR(5)),num2str(ERROR(6))});
%X = reordercats(X,{num2str(ERROR(1)),num2str(ERROR(2)),num2str(ERROR(3)),num2str(ERROR(4)),num2str(ERROR(5)),num2str(ERROR(6))});
X = categorical({num2str(ERROR(1)),num2str(ERROR(2)),num2str(ERROR(3)),num2str(ERROR(4)),num2str(ERROR(5))});
X = reordercats(X,{num2str(ERROR(1)),num2str(ERROR(2)),num2str(ERROR(3)),num2str(ERROR(4)),num2str(ERROR(5))});
A=[mean(RRR_G(:,:,1));mean(RRR_G_e(:,:,1));mean(RRR_sdma(:,:,1));mean(RRR_sdma_e(:,:,1));...
    mean(RRR_noma(:,:,1));mean(RRR_noma_e(:,:,1));mean(RRR_zf_mi(:,:,1));mean(RRR_mrt1(:,:,1))]';
bar(X,A)
legend('RS',' RS (\sigma_e=0)','SDMA','SDMA (\sigma_e=0)','NOMA','NOMA (\sigma_e=0)',...
    'multicasting','OMA-MRT')
%legend('RS','SDMA','NOMA','multicasting','OMA-MRT')
xlabel('\sigma_e^2')
ylabel('rate (bits/s/Hz)')
title({['Nt=',num2str(Nt),', M=',num2str(K),', SNR=',num2str(ERROR),', sample=',num2str(sample)], ['']})




figure(5)
hold on
X = categorical({num2str(ERROR(1)),num2str(ERROR(2)),num2str(ERROR(3)),num2str(ERROR(4)),num2str(ERROR(5))});
X = reordercats(X,{num2str(ERROR(1)),num2str(ERROR(2)),num2str(ERROR(3)),num2str(ERROR(4)),num2str(ERROR(5))});
A=[mean(RRR_G(:,:,1));mean(RRR_G_e(:,:,1));mean(RRR_zf(:,:,1));mean(RRR_zf_e(:,:,1));...
    mean(RRR_mrt(:,:,1));mean(RRR_mrt_e(:,:,1));mean(RRR_zf_si(:,:,1));mean(RRR_zf_sp(:,:,1));mean(RRR_mrt1(:,:,1))]';
bar(X,A)
legend('RS',' RS (\sigma_e=0)','RS-ZF','RS-ZF (\sigma_e=0)',...
    'RS-MRT','RS-MRT (\sigma_e=0)','SDMA-ZF','SDMA-ZF (\sigma_e=0)','OMA-MRT')

%legend('RS','RS-ZF','RS-MRT','SDMA-ZF','OMA-MRT')
xlabel('\sigma_e^2')
ylabel('rate (bits/s/Hz)')
title({['Nt=',num2str(Nt),', M=',num2str(K),', SNR=',num2str(ERROR),', sample=',num2str(sample)], ['']})




