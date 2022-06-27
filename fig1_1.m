clear all


Nt=2;%nuber of transmitter antenna
K=2;%number of user
M=1;
SNR=[0:6:30];
%Pt=10^(SNR/10);%transmit signal power

sigma_e1=0.2;
sigma_e=randn(K,1)*0+0;%M X 1 vector
result_set_e=zeros(40);
result_nc_set_e=zeros(40);
result_zf_set_e=zeros(40);
%gam=0;
beta=(sigma_e1)^2;

sample=20;
for i=1:sample
    
    i
    
    H_h=sqrt(1-beta)*(randn(Nt,K)+1i*randn(Nt,K))/sqrt(2);%Rayleigh h1
    if norm(H_h(:,1))>=norm(H_h(:,2))
        h1(:,i)=H_h(:,1); h2(:,i)=H_h(:,2);
    else
        h1(:,i)=H_h(:,2); h2(:,i)=H_h(:,1);
    end
    H_h(:,1)=h1(:,i); H_h(:,2)=h2(:,i);
    
    
    rho=1-abs(h1(:,i)'/norm(h1(:,i))*h2(:,i)/norm(h2(:,i)))^2;
    %     Rho(i)=rho;
    
    for j=1:1:length(SNR)
        
        Pt=10^(SNR(j)/10)
        N0=1+beta*Pt;%gaussina noise
        N02=1;%without channel error
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Optimization GMI
        %%%%%%%%%%%Rate splitting%%%%%%%%%%%%%%
        [P_G,pc_G,pk_G,GMI(i,j),SR(i,j)]=GMI_RS(Nt,K,H_h,H_h,Pt,M,N0,sigma_e*0);%proposed
        
        %%%%%%%%%%%%no_common_part(=SDMA)%%%%%%%%%%%%%%
        [P_sdma,pc_sdma,pk_sdma,GMI_sdma(i,j),SR_sdma(i,j)]=GMI_SDMA(Nt,K,H_h,H_h,Pt,M,N0,sigma_e*0);
        
        %%%%%%%%%Zero_Forcing%%%%%%%%%%
        [P_zf,pc_zf,pk_zf,GMI_zf(i,j),SR_zf(i,j)]=GMI_RS_ZF(Nt,K,H_h,H_h,Pt,M,N0,sigma_e*0);
        
        %%%%%%%%%MRT%%%%%%%%%%
        [P_mrt,pc_mrt,pk_mrt,GMI_mrt(i,j),SR_mrt(i,j)]=GMI_RS_MRT(Nt,K,H_h,H_h,Pt,M,N0,sigma_e*0);
        
        %%%%%%NOMA_two user
        [P_noma,pc_noma,pk_noma,GMI_noma(i,j),SR_noma(i,j)]=GMI_NOMA_two_user(Nt,K,H_h,H_h,Pt,M,N0,sigma_e*0);
        
        %================================no information about channel error============================
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Optimization GMI
        
        %%%%%%%%%%%%Rate splitting%%%%%%%%%%%%%%
        [P_G_e,pc_G_e,pk_G_e,GMI_e(i,j),SR_x(i,j)]=GMI_RS(Nt,K,H_h,H_h,Pt,M,N02,sigma_e*0);%proposed
        
        %%%%%%%%%%%no_common_part(=SDMA)%%%%%%%%%%%%%%
        [P_sdma_e,pc_sdma_e,pk_sdma_e,GMI_sdma_x(i,j),SR_sdma_x(i,j)]=GMI_SDMA(Nt,K,H_h,H_h,Pt,M,N02,sigma_e*0);
        
        %%%%%%%%%Zero_Forcing%%%%%%%%%%
        [P_zf_e,pc_zf_e,pk_zf_e,GMI_zf_x(i,j),SR_zf_x(i,j)]=GMI_RS_ZF(Nt,K,H_h,H_h,Pt,M,N02,sigma_e*0);
        
        %%%%%%%%MRT%%%%%%%%%%
        [P_mrt_e,pc_mrt_e,pk_mrt_e,GMI_mrt_x(i,j),SR_mrt_x(i,j)]=GMI_RS_MRT(Nt,K,H_h,H_h,Pt,M,N02,sigma_e*0);
        
        %%%NOMA
        [P_noma_e,pc_noma_e,pk_noma_e,GMI_noma_x(i,j),SR_noma_x(i,j)]=GMI_NOMA_two_user(Nt,K,H_h,H_h,Pt,M,N02,sigma_e*0);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%Caculation achevable Rate
        [A,B,C,D]=cal_ABCD(Nt,K,Pt,N0,H_h,sigma_e*0);
        %%%%%%%%%%%%%%Rate splitting%%%%%%%%%%%%%%%%
        SR_e(i,j)=cal_GMI(K,A,B,C,D,P_G_e);
        %%%%%%%%%%%%%%no_common_part(=SDMA)%%%%%%%%%%%%%%
        SR_sdma_e(i,j)=cal_GMI(K,A,B,C,D,P_sdma_e);
        %%%%%%%%%Zero_Forcing%%%%%%%%%%
        SR_zf_e(i,j)=cal_GMI(K,A,B,C,D,P_zf_e);
        %%%%%%%%MRT%%%%%%%%%%
        SR_mrt_e(i,j)=cal_GMI(K,A,B,C,D,P_mrt_e);
        %%%%%%%%NOMA%%%%%%%%%%
        SR_noma_e(i,j)=cal_GMI(K,A,B,C,D,P_noma_e);
        
        
        
        %%%%Zero forcing
        %%%%%%%%%RS_NOMA
        [MA_Ri(i,j),t_Ri(i,j), P1_Ri(i,j),P2_Ri(i,j), Pc_Ri(i,j),Rs_Ri(i,j),fc_Ri(:,i,j),t_Ni(i,j),Rs_Ni(i,j),fc_Ni(:,i,j)]=RSNOMA_sdr_imperfect(Pt,h1(:,i),h2(:,i),rho,beta);
        [MA_Rp(i,j),t_Rp(i,j), P1_Rp(i,j),P2_Rp(i,j), Pc_Rp(i,j),Rs_Rp(i,j),fc_Rp(:,i,j),t_Np(i,j),Rs_Np(i,j),fc_Np(:,i,j)]=RSNOMA_sdr_imperfect(Pt,h1(:,i),h2(:,i),rho,0);
        %%%%%%%%%%%%%%%%%%%%%%%%%SDMA
        [MA_Si(i,j),t_Si(i,j), P1_Si(i,j),P2_Si(i,j), Pc_Si(i,j),Rs_Si(i,j),fc_Si(:,i,j)]=SDMA_imperfect(Pt,h1(:,i),h2(:,i),rho,beta);
        [MA_Sp(i,j),t_Sp(i,j), P1_Sp(i,j),P2_Sp(i,j), Pc_Sp(i,j),Rs_Sp(i,j),fc_Sp(:,i,j)]=SDMA_perfect(Pt,h1(:,i),h2(:,i),rho);
        %%%%%%%MRT
        [MA_Oim1(i,j),t_Oim1(i,j), P1_Oim1(i,j),P2_Oim1(i,j), Pc_Oim1(i,j),Rs_Oim1(i,j),fc_Oim1(:,i,j)]=OMA_MRT_one_imperfect(Pt,h1(:,i),h2(:,i),rho,beta);
        [MA_Oim(i,j),t_Oim(i,j), P1_Oim(i,j),P2_Oim(i,j), Pc_Oim(i,j),Rs_Oim(i,j),fc_Oim(:,i,j)]=OMA_MRT_imperfect(Pt,h1(:,i),h2(:,i),rho,beta);
        %%%%%%multicasting
        [MA_mi(i,j),t_mi(i,j), P1_mi(i,j),P2_mi(i,j), Pc_mi(i,j),Rs_mi(i,j),fc_mi(:,i,j)]=multi_sdr_imperfect(Pt,h1(:,i),h2(:,i),rho,beta);
        [MA_mp(i,j),t_mp(i,j), P1_mp(i,j),P2_mp(i,j), Pc_mp(i,j),Rs_mp(i,j),fc_mp(:,i,j)]=multi_sdr_perfect(Pt,h1(:,i),h2(:,i),rho);
        
        %%%%%%%%%%%%%%%%%%=====caculate GMI=======
        %%%%%RS NOMA
        Rs__Ri(i,j)=cal_GMI_1(Pt,t_Ri(i,j), P1_Ri(i,j),P2_Ri(i,j), Pc_Ri(i,j),fc_Ri(:,i,j),h1(:,i),h2(:,i),rho,beta);
        Rs__Rp(i,j)=cal_GMI_1(Pt,t_Rp(i,j), P1_Rp(i,j),P2_Rp(i,j), Pc_Rp(i,j),fc_Rp(:,i,j),h1(:,i),h2(:,i),rho,beta);
        Rs__Ni(i,j)=cal_GMI_1(Pt,t_Ni(i,j), t_Ni(i,j)*Pt,0, (1-t_Ni(i,j))*Pt,fc_Ni(:,i,j),h1(:,i),h2(:,i),rho,beta);
        Rs__Np(i,j)=cal_GMI_1(Pt,t_Np(i,j), t_Np(i,j)*Pt,0, (1-t_Np(i,j))*Pt,fc_Np(:,i,j),h1(:,i),h2(:,i),rho,beta);
        %%%%%%%%%%%%%%%%%%%%%%%%%SDMA
        Rs__Si(i,j)=cal_GMI_1(Pt,t_Si(i,j), P1_Si(i,j),P2_Si(i,j), Pc_Si(i,j),fc_Si(:,i,j),h1(:,i),h2(:,i),rho,beta);
        Rs__Sp(i,j)=cal_GMI_1(Pt,t_Sp(i,j), P1_Sp(i,j),P2_Sp(i,j), Pc_Sp(i,j),fc_Sp(:,i,j),h1(:,i),h2(:,i),rho,beta);
        %%%%%%%%%%%%%%%multi
        Rs__mi(i,j)=cal_GMI_1(Pt,t_mi(i,j), P1_mi(i,j),P2_mi(i,j), Pc_mi(i,j),fc_mi(:,i,j),h1(:,i),h2(:,i),rho,beta);
        Rs__mp(i,j)=cal_GMI_1(Pt,t_mp(i,j), P1_mp(i,j),P2_mp(i,j), Pc_mp(i,j),fc_mp(:,i,j),h1(:,i),h2(:,i),rho,beta);
        
        

        %%MMSE


        [pc_1,pk_1]=Alg_1(Nt,K,H_h,H_h,Pt,M,N0);
        ESR(i,j)=Alg_1_cal_ESR(pc_1,pk_1,H_h,0,M,K,N0);%paper        
        
        [pc_con,pk_con,Bound(i,j)]=Alg_1_con(Nt,K,H_h,Pt,N0,sigma_e(1)^2*0);
        ESR_con(i,j)=Alg_1_cal_ESR(pc_con,pk_con,H_h,0,M,K,N0);%paper
        
        [pc_noRS,pk_noRS]=Alg_1_noRS2(Nt,K,H_h,H_h,Pt,M,N0);
        ESR_noRS(i,j)=Alg_1_cal_ESR(pc_noRS,pk_noRS,H_h,0,M,K,N0);%paper
        
    end
    
end
figure(1)
hold on
%ahevable rate
%RS imperfect with error info&%imperfect with no error info
plot(SNR,mean(SR),'-sr','LineWidth',2,'MarkerSize',8);plot(SNR,mean(SR_e),'--r>','LineWidth',1.5,'MarkerSize',8)

plot(SNR,mean(SR_sdma),'b-.s','LineWidth',1.5,'MarkerSize',8); plot(SNR,mean(SR_sdma_e),'b--^','LineWidth',1.5,'MarkerSize',8)

plot(SNR,mean(SR_zf),'k:s','LineWidth',1.5,'MarkerSize',7);
plot(SNR,mean(SR_zf_e),'k-.<','LineWidth',1.5,'MarkerSize',6);

plot(SNR,mean(SR_mrt),'y-d','LineWidth',1.5,'MarkerSize',8);
plot(SNR,mean(SR_mrt_e),'y--<','LineWidth',1.5,'MarkerSize',8);

plot(SNR,mean(SR_noma),'m--d','LineWidth',1.5,'MarkerSize',8)
plot(SNR,mean(SR_noma_e),'m--<','LineWidth',1.5,'MarkerSize',8)



%%%%ZeroForcing
%%RS
plot(SNR,mean(Rs__Ri),'k','LineWidth',1.5,'MarkerSize',4)
plot(SNR,mean(Rs__Rp),'k--o','LineWidth',1.5,'MarkerSize',8)
%%NOMA
plot(SNR,mean(Rs__Ni),'g','LineWidth',1.5,'MarkerSize',8)
plot(SNR,mean(Rs__Np),'g--o','LineWidth',1.5,'MarkerSize',8)
%SDMA
plot(SNR,mean(Rs__Si),'b','LineWidth',1.5,'MarkerSize',8)
plot(SNR,mean(Rs__Sp),'b--o','LineWidth',1.5,'MarkerSize',8)

%MRT
plot(SNR,mean(Rs_Oim1),'c--*','LineWidth',1.5,'MarkerSize',8)
plot(SNR,mean(Rs_Oim),'c--+','LineWidth',1.5,'MarkerSize',8)

%multicasting
plot(SNR,mean(Rs__mi),'m');
plot(SNR,mean(Rs__mp),'m--o');

legend('RS',' RS_{no error info}','SDMA','SDMA_{no error info}','RS-ZF','RS-ZF_{no error info}',...
    'RS-MRT','RS-MRT_{no error info}','NOMA','NOMA_{no error info}',...
    'RS-ZF','RS-ZF_{no error info}','NOMA-ZF','NOMA-ZF_{no error info}','SDMA-ZF','SDMA-ZF_{no error info}',...
    'OMA-MRT1','OMA-MRT2','multicasting','multicasting_{no error info}')

xlabel('SNR')
ylabel('rate (bits/s/Hz)')
title(['Nt=',num2str(Nt),', M=',num2str(K),', \sigma_{e}=',num2str(sigma_e1(1)),', sample=',num2str(sample)])

hold on 
%MMSE
plot(SNR,mean(ESR),'m')
plot(SNR,mean(ESR_con),'m--o')

