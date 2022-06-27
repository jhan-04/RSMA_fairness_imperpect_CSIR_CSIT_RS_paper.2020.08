%clear all


Nt=2;%nuber of transmitter antenna
K=2;%number of user
M=1;

SNR=[5:6:35];
%Pt=10^(SNR/10);%transmit signal power

sigma_e1=sqrt(0.3);
beta=(sigma_e1)^2;

sigma_e=randn(K,1)*0+0;%M X 1 vector
result_set_e=zeros(40);
result_nc_set_e=zeros(40);
result_zf_set_e=zeros(40);
%gam=0;

sample=100;
for i=1:sample
    
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
    %     %     Rho(i)=rho;
    
    for j=1:1:length(SNR)
        
        
        snr=SNR(j)
        Pt=10^(snr/10);
        
        
        
        
        rho=1-abs(h1(:,i)'/norm(h1(:,i))*h2(:,i)/norm(h2(:,i)))^2;
        N0=1+beta*Pt;%gaussina noise
        N02=1;%without channel error
        
        
        %%%%%%%%%%%Rate splitting%%%%%%%%%%%%%%
        [P_G,pc_G,pk_G,GMI(i,j),SR(i,j)]=GMI_RS(Nt,K,H_h,H_h,Pt,M,N0,sigma_e*0);%proposed
        
        %%%%%%%%%Zero_Forcing%%%%%%%%%%
        [P_zf,pc_zf,pk_zf,GMI_zf(i,j),SR_zf(i,j)]=GMI_RS_ZF(Nt,K,H_h,H_h,Pt,M,N0,sigma_e*0);
        
        %%%%%%%%%MRT%%%%%%%%%%
        [P_mrt,pc_mrt,pk_mrt,GMI_mrt(i,j),SR_mrt(i,j)]=GMI_RS_MRT(Nt,K,H_h,H_h,Pt,M,N0,sigma_e*0);
        
 
        %================================no information about channel error============================
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Optimization GMI
        
        %%%%%%%%%%%%Rate splitting%%%%%%%%%%%%%%
        [P_G_e,pc_G_e,pk_G_e,GMI_e(i,j),SR_x(i,j)]=GMI_RS(Nt,K,H_h,H_h,Pt,M,N02,sigma_e*0);%proposed
        %%%%%%%%%Zero_Forcing%%%%%%%%%%
        [P_zf_e,pc_zf_e,pk_zf_e,GMI_zf_x(i,j),SR_zf_x(i,j)]=GMI_RS_ZF(Nt,K,H_h,H_h,Pt,M,N02,sigma_e*0);
        %%%%%%%%MRT%%%%%%%%%%
        [P_mrt_e,pc_mrt_e,pk_mrt_e,GMI_mrt_x(i,j),SR_mrt_x(i,j)]=GMI_RS_MRT(Nt,K,H_h,H_h,Pt,M,N02,sigma_e*0);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%Caculation achevable Rate
        [A,B,C,D]=cal_ABCD(Nt,K,Pt,N0,H_h,sigma_e*0);
        %%%%%%%%%%%%%%Rate splitting%%%%%%%%%%%%%%%%
        SR_e(i,j)=cal_GMI(K,A,B,C,D,P_G_e);
        %%%%%%%%%Zero_Forcing%%%%%%%%%%
        SR_zf_e(i,j)=cal_GMI(K,A,B,C,D,P_zf_e);
        %%%%%%%%MRT%%%%%%%%%%
        SR_mrt_e(i,j)=cal_GMI(K,A,B,C,D,P_mrt_e);

        
    end
    
end

figure(4)
hold on


plot(SNR,mean(SR),'-sr','LineWidth',2,'MarkerSize',8)
plot(SNR,mean(SR_e),'--r>','LineWidth',1.5,'MarkerSize',8);


plot(SNR,mean(SR_zf),'k:s','LineWidth',1.5,...
    'MarkerSize',7)
plot(SNR,mean(SR_zf_e),'k-.<','LineWidth',1.5,...
    'MarkerSize',6)

plot(SNR,mean(SR_mrt),'g-d','LineWidth',1.5,...
    'MarkerSize',8)
plot(SNR,mean(SR_mrt_e),'g--<','LineWidth',1.5,...
    'MarkerSize',8)




legend('RS',' RS (\sigma_e=0)','RS-ZF',' RS-ZF (\sigma_e=0)','RS-MRT',' RS-MRT (\sigma_e=0)')
xlabel('SNR (dB)')
ylabel('Sum-rate (bits/s/Hz)')
title({['Nt=',num2str(Nt),', M=',num2str(K),', \sigma_e=',num2str(sigma_e1),', sample=',num2str(sample)], ['']})



hold on
plot(SNR,mean(RRR_G(:,:,1)),'-sr','LineWidth',2,'MarkerSize',8);
plot(SNR,mean(RRR_zf(:,:,1)),'k:s','LineWidth',1.5,'MarkerSize',7);
plot(SNR,mean(RRR_mrt(:,:,1)),'g-d','LineWidth',1.5,'MarkerSize',8);
legend('RS','RS-ZF','RS-MRT')

