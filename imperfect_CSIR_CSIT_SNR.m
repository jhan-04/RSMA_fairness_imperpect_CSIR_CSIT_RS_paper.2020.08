clear all



Nt=2;%nuber of transmitter antenna
K=2;%number of user
M=1;
SNR=[0:6:30];
%Pt=10^(SNR/10);%transmit signal power

sigma_e1=0.3;
sigma_e=randn(K,1)*0+0;%M X 1 vector
result_set_e=zeros(40);
result_nc_set_e=zeros(40);
result_zf_set_e=zeros(40);
%gam=0;
sigma_e2=(sigma_e1)^2;

sample=10;
for i=1:sample
    
    i
    
    H_h=sqrt(1-sigma_e1^2)*(randn(Nt,K)+1i*randn(Nt,K))/sqrt(2);%Rayleigh h1
    if norm(H_h(:,1))>=norm(H_h(:,2))
        h1(:,i)=H_h(:,1); h2(:,i)=H_h(:,2);
    else
        h1(:,i)=H_h(:,2); h2(:,i)=H_h(:,1);
    end
    H_h(:,1)=h1(:,i); H_h(:,2)=h2(:,i);
    
    
    rho=1-abs(h1(:,i)'/norm(h1(:,i))*h2(:,i)/norm(h2(:,i)))^2;
    %     Rho(i)=rho;
    
    for j=1:1:length(SNR)
        
        %%hK
        
        
        Pt=10^(SNR(j)/10);
        
        N0=1+sigma_e2*Pt;%gaussina noise
        N02=1;%without channel error
        
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Optimization GMI
%         %%%%%%%%%%%%Rate splitting%%%%%%%%%%%%%%
%         [P_G,pc_G,pk_G,GMI(i,j),SR(i,j)]=GMI_RS(Nt,K,H_h,H_h,Pt,M,N0,sigma_e*0);%proposed
%         
%         %%%%%%%%%%%%no_common_part(=SDMA)%%%%%%%%%%%%%%
%         [P_sdma,pc_sdma,pk_sdma,GMI_sdma(i,j),SR_sdma(i,j)]=GMI_SDMA(Nt,K,H_h,H_h,Pt,M,N0,sigma_e*0);
%         
%         %%%%%%%%%Zero_Forcing%%%%%%%%%%
%         [P_zf,pc_zf,pk_zf,GMI_zf(i,j),SR_zf(i,j)]=GMI_RS_ZF(Nt,K,H_h,H_h,Pt,M,N0,sigma_e*0);
%         
%         %%%%%%%%%MRT%%%%%%%%%%
%         [P_mrt,pc_mrt,pk_mrt,GMI_mrt(i,j),SR_mrt(i,j)]=GMI_RS_MRT(Nt,K,H_h,H_h,Pt,M,N0,sigma_e*0);
%         
%         
%         
%         %================================no information about channel error============================
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Optimization GMI
%         
%         %%%%%%%%%%%%Rate splitting%%%%%%%%%%%%%%
%         [P_G_e,pc_G_e,pk_G_e,GMI_e(i,j),SR_e(i,j)]=GMI_RS(Nt,K,H_h,H_h,Pt,M,N02,sigma_e*0);%proposed
% 
%         %%%%%%%%%%%no_common_part(=SDMA)%%%%%%%%%%%%%%
%        [P_sdma_e,pc_sdma_e,pk_sdma_e,GMI_sdma_e(i,j),SR_sdma_e(i,j)]=GMI_SDMA(Nt,K,H_h,H_h,Pt,M,N02,sigma_e*0);
%         
%         %%%%%%%%%Zero_Forcing%%%%%%%%%%
%         [P_zf_e,pc_zf_e,pk_zf_e,GMI_zf_e(i,j),SR_zf_e(i,j)]=GMI_RS_ZF(Nt,K,H_h,H_h,Pt,M,N02,sigma_e*0); 
%         
%         %%%%%%%%MRT%%%%%%%%%%
%        [P_mrt_e,pc_mrt_e,pk_mrt_e,GMI_mrt_e(i,j),SR_mrt_e(i,j)]=GMI_RS_MRT(Nt,K,H_h,H_h,Pt,M,N02,sigma_e*0); 
%         
%         %%%%%%%%%%%%%%%%%%%%%%%%%%Caculation achevable Rate 
%         [A,B,C,D]=cal_ABCD(Nt,K,Pt,N0,H_h,sigma_e*0);
%         %%%%%%%%%%%%%%Rate splitting%%%%%%%%%%%%%%%%
%         Rs_e(i,j)=cal_GMI(K,A,B,C,D,P_G_e);
%         %%%%%%%%%%%%%%no_common_part(=SDMA)%%%%%%%%%%%%%%
%         Rs_nc_e(i,j)=cal_GMI(K,A,B,C,D,P_sdma_e);
%         %%%%%%%%%Zero_Forcing%%%%%%%%%%
%         Rs_zf_e(i,j)=cal_GMI(K,A,B,C,D,P_zf_e);
%         %%%%%%%%MRT%%%%%%%%%%
%         Rs_mrt_e(i,j)=cal_GMI(K,A,B,C,D,P_mrt_e);
        
    
        
        
        
        %%%%%%%%%%RS_NOMA
        %%%%%%%%%RS_NOMA
        [MA_Ri(i,j),t_Ri(i,j), P1_Ri(i,j),P2_Ri(i,j), Pc_Ri(i,j),Rs_Ri(i,j),fc_Ri(:,i,j),t_Ni(i,j),Rs_Ni(i,j),fc_Ni(:,i,j)]=RSNOMA_sdr_imperfect(Pt,h1(:,i),h2(:,i),rho,sigma_e2);
        [MA_Rp(i,j),t_Rp(i,j), P1_Rp(i,j),P2_Rp(i,j), Pc_Rp(i,j),Rs_Rp(i,j),fc_Rp(:,i,j),t_Np(i,j),Rs_Np(i,j),fc_Np(:,i,j)]=RSNOMA_sdr_imperfect(Pt,h1(:,i),h2(:,i),rho,0);


        
        %%%%%%%%%%%%%%%%%%%%%%%%%SDMA
        [MA_Si(i,j),t_Si(i,j), P1_Si(i,j),P2_Si(i,j), Pc_Si(i,j),Rs_Si(i,j),fc_Si(:,i,j)]=SDMA_imperfect(Pt,h1(:,i),h2(:,i),rho,sigma_e2);
        [MA_Sp(i,j),t_Sp(i,j), P1_Sp(i,j),P2_Sp(i,j), Pc_Sp(i,j),Rs_Sp(i,j),fc_Sp(:,i,j)]=SDMA_perfect(Pt,h1(:,i),h2(:,i),rho);
        
        %%%%%%%OMA
        [MA_Oi(i,j),t_Oi(i,j), P1_Oi(i,j),P2_Oi(i,j), Pc_Oi(i,j),Rs_Oi(i,j),fc_Oi(:,i,j)]=OMA_imperfect(Pt,h1(:,i),h2(:,i),rho,sigma_e2);
        [MA_Op(i,j),t_Op(i,j), P1_Op(i,j),P2_Op(i,j), Pc_Op(i,j),Rs_Op(i,j),fc_Op(:,i,j)]=OMA_perfect(Pt,h1(:,i),h2(:,i),rho);
        %%%%%%%MRT
        
        [MA_Oim1(i,j),t_Oim1(i,j), P1_Oim1(i,j),P2_Oim1(i,j), Pc_Oim1(i,j),Rs_Oim1(i,j),fc_Oim1(:,i,j)]=OMA_MRT_one_imperfect(Pt,h1(:,i),h2(:,i),rho,sigma_e2);Rs_Oim1(i,j)
        [MA_Oim(i,j),t_Oim(i,j), P1_Oim(i,j),P2_Oim(i,j), Pc_Oim(i,j),Rs_Oim(i,j),fc_Oim(:,i,j)]=OMA_MRT_imperfect(Pt,h1(:,i),h2(:,i),rho,sigma_e2);
        
        %%%%%%multicasting
        [MA_mi(i,j),t_mi(i,j), P1_mi(i,j),P2_mi(i,j), Pc_mi(i,j),Rs_mi(i,j),fc_mi(:,i,j)]=multi_sdr_imperfect(Pt,h1(:,i),h2(:,i),rho,sigma_e2);
        [MA_mp(i,j),t_mp(i,j), P1_mp(i,j),P2_mp(i,j), Pc_mp(i,j),Rs_mp(i,j),fc_mp(:,i,j)]=multi_sdr_perfect(Pt,h1(:,i),h2(:,i),rho);
        %%KKT
        [MA_mik(i,j),t_mik(i,j), P1_mik(i,j),P2_mik(i,j), Pc_mik(i,j),Rs_mik(i,j),fc_mik(:,i,j)]=multi_kkt_imperfect(Pt,h1(:,i),h2(:,i),rho,sigma_e2);
        [MA_mpk(i,j),t_mpk(i,j), P1_mpk(i,j),P2_mpk(i,j), Pc_mpk(i,j),Rs_mpk(i,j),fc_mpk(:,i,j)]=multi_kkt_perfect(Pt,h1(:,i),h2(:,i),rho);
        
        
        %%%%%%%%%%%%%%%%%%caculate GMI
        
        Rs__Ri(i,j)=cal_GMI_1(Pt,t_Ri(i,j), P1_Ri(i,j),P2_Ri(i,j), Pc_Ri(i,j),fc_Ri(:,i,j),h1(:,i),h2(:,i),rho,sigma_e2);
        Rs__Rp(i,j)=cal_GMI_1(Pt,t_Rp(i,j), P1_Rp(i,j),P2_Rp(i,j), Pc_Rp(i,j),fc_Rp(:,i,j),h1(:,i),h2(:,i),rho,sigma_e2);
        Rs__Ni(i,j)=cal_GMI_1(Pt,t_Ni(i,j), t_Ni(i,j)*Pt,0, (1-t_Ni(i,j))*Pt,fc_Ni(:,i,j),h1(:,i),h2(:,i),rho,sigma_e2);
        Rs__Np(i,j)=cal_GMI_1(Pt,t_Np(i,j), t_Np(i,j)*Pt,0, (1-t_Np(i,j))*Pt,fc_Np(:,i,j),h1(:,i),h2(:,i),rho,sigma_e2);
        
        %%%KKT
        %         Rs__Rik(i,j)=cal_GMI_1(Pt,t_Rik(i,j), P1_Rik(i,j),P2_Rik(i,j), Pc_Rik(i,j),fc_Rik(:,i,j),h1(:,i),h2(:,i),rho,beta);
        %         Rs__Rpk(i,j)=cal_GMI_1(Pt,t_Rpk(i,j), P1_Rpk(i,j),P2_Rpk(i,j), Pc_Rpk(i,j),fc_Rpk(:,i,j),h1(:,i),h2(:,i),rho,beta);
        %         Rs__Nik(i,j)=cal_GMI_1(Pt,t_Nik(i,j), t_Nik(i,j)*Pt, 0, (1-t_Nik(i,j))*Pt,fc_Nik(:,i,j),h1(:,i),h2(:,i),rho,beta);
        %         Rs__Npk(i,j)=cal_GMI_1(Pt,t_Npk(i,j), t_Npk(i,j)*Pt, 0, (1-t_Npk(i,j))*Pt,fc_Npk(:,i,j),h1(:,i),h2(:,i),rho,beta);
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%SDMA
        Rs__Si(i,j)=cal_GMI_1(Pt,t_Si(i,j), P1_Si(i,j),P2_Si(i,j), Pc_Si(i,j),fc_Si(:,i,j),h1(:,i),h2(:,i),rho,sigma_e2);
        Rs__Sp(i,j)=cal_GMI_1(Pt,t_Sp(i,j), P1_Sp(i,j),P2_Sp(i,j), Pc_Sp(i,j),fc_Sp(:,i,j),h1(:,i),h2(:,i),rho,sigma_e2);
        
        %%%%%%%%%%%%%%%%%%%OMA
        Rs__Oi(i,j)=cal_GMI_1(Pt,t_Oi(i,j), P1_Oi(i,j),P2_Oi(i,j), Pc_Oi(i,j),fc_Oi(:,i,j),h1(:,i),h2(:,i),rho,sigma_e2);
        Rs__Op(i,j)=cal_GMI_1(Pt,t_Op(i,j), P1_Op(i,j),P2_Op(i,j), Pc_Op(i,j),fc_Op(:,i,j),h1(:,i),h2(:,i),rho,sigma_e2);
        
        %%%%%%%%%%%%%%%multi
        Rs__mi(i,j)=cal_GMI_1(Pt,t_mi(i,j), P1_mi(i,j),P2_mi(i,j), Pc_mi(i,j),fc_mi(:,i,j),h1(:,i),h2(:,i),rho,sigma_e2);
        Rs__mp(i,j)=cal_GMI_1(Pt,t_mp(i,j), P1_mp(i,j),P2_mp(i,j), Pc_mp(i,j),fc_mp(:,i,j),h1(:,i),h2(:,i),rho,sigma_e2);
        %%KKT
        Rs__mik(i,j)=cal_GMI_1(Pt,t_mik(i,j), P1_mik(i,j),P2_mik(i,j), Pc_mik(i,j),fc_mik(:,i,j),h1(:,i),h2(:,i),rho,sigma_e2);
        Rs__mpk(i,j)=cal_GMI_1(Pt,t_mpk(i,j), P1_mpk(i,j),P2_mpk(i,j), Pc_mpk(i,j),fc_mpk(:,i,j),h1(:,i),h2(:,i),rho,sigma_e2);
        
        
        
        
        
        
        fprintf('=======sample=%d===================SNR=%d==================\n',i,SNR(j))
        
        fprintf('RS_SDR_im: MA=%d , t=%d, P1=%d, P2=%d Pc=%d Rs=%d GMI=%d \n ',MA_Ri(i,j),t_Ri(i,j), P1_Ri(i,j),P2_Ri(i,j), Pc_Ri(i,j),Rs_Ri(i,j),Rs__Ri(i,j))
        fprintf('NOMA_SDR_im: t=%d, Rs=%d Rs2=%d\n\n',t_Ni(i,j),Rs_Ni(i,j),Rs__Ni(i,j))
        fprintf('RS_SDR: MA=%d , t=%d, P1=%d, P2=%d Pc=%d Rs=%d GMI=%d\n ',MA_Rp(i,j),t_Rp(i,j), P1_Rp(i,j),P2_Rp(i,j), Pc_Rp(i,j),Rs_Rp(i,j),Rs__Rp(i,j))
        fprintf('NOMA_SDR: t=%d, Rs=%d GMI=%d\n\n',t_Np(i,j),Rs_Np(i,j),Rs__Np(i,j))
        %%%
        %         fprintf('RS_KKT_im: MA=%d , t=%d, P1=%d, P2=%d Pc=%d Rs=%d GMI=%d \n ',MA_Rik(i,j),t_Rik(i,j), P1_Rik(i,j),P2_Rik(i,j), Pc_Rik(i,j),Rs_Rik(i,j),Rs__Rik(i,j))
        %         fprintf('NOMA_KKT_im: t=%d, Rs=%d Rs2=%d\n\n',t_Nik(i,j),Rs_Nik(i,j),Rs__Nik(i,j))
        %         fprintf('RS_KKT: MA=%d , t=%d, P1=%d, P2=%d Pc=%d Rs=%d GMI=%d\n ',MA_Rpk(i,j),t_Rpk(i,j), P1_Rpk(i,j),P2_Rpk(i,j), Pc_Rpk(i,j),Rs_Rpk(i,j),Rs__Rpk(i,j))
        %         fprintf('NOMA_KKT: t=%d, Rs=%d GMI=%d\n\n',t_Npk(i,j),Rs_Npk(i,j),Rs__Npk(i,j))
        
        fprintf('SDMA_im: MA=%d , t=%d, P1=%d, P2=%d Pc=%d Rs=%d GMI=%d\n ',MA_Si(i,j),t_Si(i,j), P1_Si(i,j),P2_Si(i,j), Pc_Si(i,j),Rs_Si(i,j),Rs__Si(i,j))
        fprintf('SDMA: MA=%d , t=%d, P1=%d, P2=%d Pc=%d Rs=%d GMI=%d\n ',MA_Sp(i,j),t_Sp(i,j), P1_Sp(i,j),P2_Sp(i,j), Pc_Sp(i,j),Rs_Sp(i,j),Rs__Sp(i,j))
        
        fprintf('OMA_im: MA=%d , t=%d, P1=%d, P2=%d Pc=%d Rs=%d GMI=%d \n ',MA_Oi(i,j),t_Oi(i,j), P1_Oi(i,j),P2_Oi(i,j), Pc_Oi(i,j),Rs_Oi(i,j),Rs__Oi(i,j))
        fprintf('OMA: MA=%d , t=%d, P1=%d, P2=%d Pc=%d Rs=%d GMI=%d\n ',MA_Op(i,j),t_Op(i,j), P1_Op(i,j),P2_Op(i,j), Pc_Op(i,j),Rs_Op(i,j),Rs__Op(i,j))
        
        fprintf('multi_SDR_im: MA=%d , t=%d, P1=%d, P2=%d Pc=%d Rs=%d GMI=%d\n ',MA_mi(i,j),t_mi(i,j), P1_mi(i,j),P2_mi(i,j), Pc_mi(i,j),Rs_mi(i,j),Rs__mi(i,j))
        fprintf('multi_SDR: MA=%d , t=%d, P1=%d, P2=%d Pc=%d Rs=%d GMI=%d\n ',MA_mp(i,j),t_mp(i,j), P1_mp(i,j),P2_mp(i,j), Pc_mp(i,j),Rs_mp(i,j),Rs__mp(i,j))
        
        fprintf('multi_KKT_im: MA=%d , t=%d, P1=%d, P2=%d Pc=%d Rs=%d GMI=%d\n ',MA_mik(i,j),t_mik(i,j), P1_mik(i,j),P2_mik(i,j), Pc_mik(i,j),Rs_mik(i,j),Rs__mik(i,j))
        fprintf('multi_KKT: MA=%d , t=%d, P1=%d, P2=%d Pc=%d Rs=%d GMI=%d\n ',MA_mpk(i,j),t_mpk(i,j), P1_mpk(i,j),P2_mpk(i,j), Pc_mpk(i,j),Rs_mpk(i,j),Rs__mpk(i,j))
        
        
        
    end
    
    
    
    
    
end

%i=1:sample;
hold on
%ahevable rate
%imperfect with error info
plot(SNR,mean(GMI_L),'-sr','LineWidth',2,...
    'MarkerSize',8)
plot(SNR,mean(GMI_L_nc),'b-.s','LineWidth',1.5,...
    'MarkerSize',8)
plot(SNR,mean(GMI_zf),'k:s','LineWidth',1.5,...
    'MarkerSize',7)
plot(SNR,mean(GMI_mrt),'y-d','LineWidth',1.5,...
    'MarkerSize',8)
%imperfect with no error info
plot(SNR,mean(ach_rate_e),'--r>','LineWidth',1.5,...
    'MarkerSize',8)
plot(SNR,mean(ach_rate_nc_e),'b--^','LineWidth',1.5,...
    'MarkerSize',8)
plot(SNR,mean(ach_rate_zf_e),'k-.<','LineWidth',1.5,...
    'MarkerSize',6)
plot(SNR,mean(ach_rate_mrt_e),'y--<','LineWidth',1.5,...
    'MarkerSize',8)

%legend('RS','SDMA','RS-ZF','RS-MRT',' RS_{no error info}','SDMA_{no error info}','RS-ZF_{no error info}','RS-MRT_{no error info}')


%%%%ZeroForcing

%%RS
plot(SNR,mean(Rs__Ri),'k','LineWidth',1.5,...
    'MarkerSize',4)
plot(SNR,mean(Rs__Rp),'k--o','LineWidth',1.5,...
    'MarkerSize',8)
%%NOMA
plot(SNR,mean(Rs__Ni),'g','LineWidth',1.5,...
    'MarkerSize',8)
plot(SNR,mean(Rs__Np),'g--o','LineWidth',1.5,...
    'MarkerSize',8)
%plot(SNR,mean(Rs_Nim),'b--d')
%SDMA
plot(SNR,mean(Rs__Si),'b','LineWidth',1.5,...
    'MarkerSize',8)
plot(SNR,mean(Rs__Sp),'b--o','LineWidth',1.5,...
    'MarkerSize',8)
%OMA
plot(SNR,mean(Rs__Oi),'c')
plot(SNR,mean(Rs__Op),'c--o','LineWidth',1.5,...
    'MarkerSize',8)
%MRT
plot(SNR,mean(Rs_Oim1),'c--*','LineWidth',1.5,...
    'MarkerSize',8)
plot(SNR,mean(Rs_Oim),'c--+','LineWidth',1.5,...
    'MarkerSize',8)

%multicasting
plot(SNR,mean(Rs__mi),'m')
plot(SNR,mean(Rs__mp),'m--o')
%legend('RS-ZF','RS-ZF_{no error info}','NOMA-ZF','NOMA-ZF_{no error info}','SDMA-ZF','SDMA-ZF_{no error info}','OMA-ZF','OMA-ZF_{no error info}','OMA-MRT1','OMA-MRT2','multicasting','multicasting_{no error info}')

legend('RS','SDMA','RS-ZF','RS-MRT',' RS_{no error info}','SDMA_{no error info}','RS-ZF_{no error info}','RS-MRT_{no error info}','RS-ZF','RS-ZF_{no error info}','NOMA-ZF','NOMA-ZF_{no error info}','SDMA-ZF','SDMA-ZF_{no error info}','OMA-ZF','OMA-ZF_{no error info}','OMA-MRT1','OMA-MRT2','multicasting','multicasting_{no error info}')

xlabel('SNR')
ylabel('rate (bits/s/Hz)')
title(['Nt=',num2str(Nt),', M=',num2str(K),', \sigma_{e}=',num2str(sigma_e1(1)),', sample=',num2str(sample)])



%pp=sqrt(Pt)*p_max/norm(p_max);

