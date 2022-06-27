% clear all


Nt=2;%nuber of transmitter antenna
K=2;%number of user
M=1;

SNR=[-5:10:35];
%Pt=10^(SNR/10);%transmit signal power

sigma_e1=sqrt(0);
beta=(sigma_e1)^2;

sigma_e=randn(K,1)*0+0;%M X 1 vector
result_set_e=zeros(40);
result_nc_set_e=zeros(40);
result_zf_set_e=zeros(40);
%gam=0;

sample=100;
for i=28:sample
    
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
%     
    for j=1:1:length(SNR)
        
        
        snr=SNR(j)
        Pt=10^(snr/10);
        
        
        
        
        rho=1-abs(h1(:,i)'/norm(h1(:,i))*h2(:,i)/norm(h2(:,i)))^2;
        N0=1+beta*Pt;%gaussina noise
        N02=1;%without channel error
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Optimization GMI
        [RRR_G(i,j,:),RRR_sdma(i,j,:),RRR_zf(i,j,:),RRR_mrt(i,j,:),RRR_noma(i,j,:)]=CAL_SR_poroposed(Nt,K,H_h,H_h,Pt,M,N0,sigma_e*0);
        %%%% ================================no information about channel error=
        [RRR_G_e(i,j,:),RRR_sdma_e(i,j,:),RRR_zf_e(i,j,:),RRR_mrt_e(i,j,:),RRR_noma_e(i,j,:)]=CAL_SR_poroposed_noinfo(Nt,K,H_h,H_h,Pt,M,N0,N02,sigma_e*0);
        %Zero forcing
        [RRR_zf_i(i,j,:),RRR_zf_p(i,j,:),RRR_zf_ni(i,j,:),RRR_zf_np(i,j,:),RRR_zf_si(i,j,:),RRR_zf_sp(i,j,:),...
            RRR_zf_mi(i,j,:),RRR_zf_mp(i,j,:),RRR_mrt1(i,j,:),RRR_mrt2(i,j,:)]=CAL_SR_ZF_poroposed_twoUser(H_h,Pt,beta);
        

        %%%%%NOMA_two user
        
        [A,B,C,D]=cal_ABCD(Nt,K,Pt,N0,H_h,sigma_e*0);
        [X_N,result_set_N]=cal_X_NOMA2(Nt,K,H_h,A,B,C,D);
        [p_N(:,i,j)]=find_p_NOMA(K,Pt,A,B,C,D,X_N,H_h);
        [RRR_noma2(i,j,:)]=cal_GMI(K,A,B,C,D,p_N(:,i,j));
        
        %%%
[P_noma,pc_noma,pk_noma,GMI_noma,SR_noma]=GMI_NOMA_two_user(Nt,K,H_h,H_h,Pt,M,N0,sigma_e*0);
[SR_noma,Rs_c_noma,Rs_k_noma]=Cal_Rate(pc_noma,pk_noma,H_h,K,N0);
RRR_noma3(i,j,:)=[SR_noma,Rs_c_noma,Rs_k_noma]';

        RRR_noma3(i,j,:)
        RRR_noma2(i,j,:)


        [AA,BB,CC,DD]=cal_ABCD(Nt,K,Pt,N02,H_h,sigma_e*0);
        [X_N_e,result_set_N]=cal_X_NOMA2(Nt,K,H_h,AA,BB,CC,DD);
        [p_N_e(:,i,j)]=find_p_NOMA(K,Pt,AA,BB,CC,DD,X_N_e,H_h);
        [RRR_noma_e2(i,j,:)]=cal_GMI(K,A,B,C,D,p_N_e(:,i,j));
        
        
[Rs_x,Rs_x1,Rs_x2]=OMA_MRT_imperfect(Pt,h1(:,i),h2(:,i),rho,beta);
RRR_mrt3(i,j,:)=[Rs_x,0,Rs_x1,Rs_x2]';
        



        RRR_noma(i,j,:);
        RRR_mrt1(i,j,:) ;
        %MMSE
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

plot(SNR,mean(RRR_mrt(:,:,1)),'y-d','LineWidth',1.5,'MarkerSize',8);plot(SNR,mean(RRR_mrt_e(:,:,1)),'y--<','LineWidth',1.5,'MarkerSize',8);

plot(SNR,mean(RRR_noma2(:,:,1)),'m--d','LineWidth',1.5,'MarkerSize',8);plot(SNR,mean(RRR_noma_e2(:,:,1)),'m--<','LineWidth',1.5,'MarkerSize',8);


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
% plot(SNR,mean(ESR),'m')
% plot(SNR,mean(ESR_con),'m--o')

figure(2)
hold on
plot(SNR,mean(RRR_G(:,:,1)),'-sr','LineWidth',2,'MarkerSize',8);plot(SNR,mean(RRR_G_e(:,:,1)),'--r>','LineWidth',1.5,'MarkerSize',8);
plot(SNR,mean(RRR_sdma(:,:,1)),'b-.s','LineWidth',1.5,'MarkerSize',8);plot(SNR,mean(RRR_sdma_e(:,:,1)),'b--^','LineWidth',1.5,'MarkerSize',8);
plot(SNR,mean(RRR_noma2(:,:,1)),'k--s','LineWidth',1.5,'MarkerSize',8);plot(SNR,mean(RRR_noma_e2(:,:,1)),'k--<','LineWidth',1.5,'MarkerSize',8);
plot(SNR,mean(RRR_zf_mi(:,:,1)),'m-o','LineWidth',1.5,'MarkerSize',8);
plot(SNR,mean(RRR_mrt1(:,:,1)),'c--*','LineWidth',1.5,'MarkerSize',8);
legend('RS',' RS (\sigma_e=0)','SDMA','SDMA (\sigma_e=0)','NOMA','NOMA (\sigma_e=0)',...
    'multicasting','OMA-MRT')
xlabel('SNR (dB)')
ylabel('rate (bits/s/Hz)')
title({['Nt=',num2str(Nt),', M=',num2str(K),', \sigma_e=',num2str(sigma_e1),', sample=',num2str(sample)], ['']})



figure(3)
hold on
plot(SNR,mean(RRR_G(:,:,1)),'-sr','LineWidth',2,'MarkerSize',8);plot(SNR,mean(RRR_G_e(:,:,1)),'--r>','LineWidth',1.5,'MarkerSize',8);
plot(SNR,mean(RRR_zf(:,:,1)),'k:s','LineWidth',1.5,'MarkerSize',7);plot(SNR,mean(RRR_zf_e(:,:,1)),'k-.<','LineWidth',1.5,'MarkerSize',6);
plot(SNR,mean(RRR_mrt(:,:,1)),'g-d','LineWidth',1.5,'MarkerSize',8);plot(SNR,mean(RRR_mrt_e(:,:,1)),'g--<','LineWidth',1.5,'MarkerSize',8);
plot(SNR,mean(RRR_zf_si(:,:,1)),'b','LineWidth',1.5,'MarkerSize',8);plot(SNR,mean(RRR_zf_sp(:,:,1)),'b--o','LineWidth',1.5,'MarkerSize',8);
plot(SNR,mean(RRR_mrt1(:,:,1)),'c--*','LineWidth',1.5,'MarkerSize',8);
legend('RS',' RS (\sigma_e=0)','RS-ZF','RS-ZF (\sigma_e=0)',...
    'RS-MRT','RS-MRT (\sigma_e=0)','SDMA-ZF','SDMA-ZF (\sigma_e=0)','OMA-MRT')
xlabel('SNR (dB)')
ylabel('rate (bits/s/Hz)')
title({['Nt=',num2str(Nt),', M=',num2str(K),', \sigma_e=',num2str(sigma_e1),', sample=',num2str(sample)], ''})


figure(4)
hold on
plot(SNR,mean(RRR_G(:,:,1)),'-sr','LineWidth',2,'MarkerSize',8);plot(SNR,mean(RRR_G_e(:,:,1)),'--r>','LineWidth',1.5,'MarkerSize',8);
plot(SNR,mean(RRR_sdma(:,:,1)),'b-.s','LineWidth',1.5,'MarkerSize',8);plot(SNR,mean(RRR_sdma_e(:,:,1)),'b--^','LineWidth',1.5,'MarkerSize',8);
plot(SNR,mean(RRR_noma2(:,:,1)),'k--s','LineWidth',1.5,'MarkerSize',8);plot(SNR,mean(RRR_noma_e2(:,:,1)),'k--<','LineWidth',1.5,'MarkerSize',8);
plot(SNR,mean(RRR_zf_mi(:,:,1)),'m-o','LineWidth',1.5,'MarkerSize',8);
plot(SNR,mean(RRR_mrt2(:,:,1)),'c--*','LineWidth',1.5,'MarkerSize',8);
plot(SNR,mean(RRR_zf_si(:,:,1)),'g--s','LineWidth',1.5,'MarkerSize',8);plot(SNR,mean(RRR_zf_sp(:,:,1)),'g--<','LineWidth',1.5,'MarkerSize',8);
legend('RS',' RS (\sigma_e=0)','SDMA','SDMA (\sigma_e=0)','NOMA','NOMA (\sigma_e=0)',...
    'Multicasting','OMA-MRT','SDMA-ZF','SDMA-ZF (\sigma_e=0)')
xlabel('SNR (dB)')
ylabel('Sum-rate (bits/s/Hz)')
title({['Nt=',num2str(Nt),', M=',num2str(K),', \sigma_e=',num2str(sigma_e1),', sample=',num2str(sample)], ['']})


% %%$calculatie rate for user
% % RRR_G
% R1_G=RRR_G(:,:,2)./2+RRR_G(:,:,3);
% R2_G=RRR_G(:,:,2)./2+RRR_G(:,:,4);
% % RRR_sdma
% R1_sdma=RRR_sdma(:,:,2)./2+RRR_sdma(:,:,3);
% R2_sdma=RRR_sdma(:,:,2)./2+RRR_sdma(:,:,4);
% % RRR_zf
% R1_zf=RRR_zf(:,:,2)./2+RRR_zf(:,:,3);
% R2_zf=RRR_zf(:,:,2)./2+RRR_zf(:,:,4);
% % RRR_mrt
% R1_mrt=RRR_mrt(:,:,2)./2+RRR_mrt(:,:,3);
% R2_mrt=RRR_mrt(:,:,2)./2+RRR_mrt(:,:,4);
% % RRR_noma
% R1_noma=RRR_noma(:,:,2)./2+RRR_noma(:,:,3);
% R2_noma=RRR_noma(:,:,2)./2+RRR_noma(:,:,4);
% % RRR_G_e
% % RRR_sdma_e
% % RRR_zf_e
% % RRR_mrt_e
% % RRR_noma_e
% %%%
% % RRR_zf_i
% % RRR_zf_p
% % RRR_zf_ni
% % RRR_zf_np
% % RRR_zf_si
% R1_zf_si=RRR_zf_si(:,:,2)./2+RRR_zf_si(:,:,3);
% R2_zf_si=RRR_zf_si(:,:,2)./2+RRR_zf_si(:,:,4);
% % RRR_zf_sp
% % RRR_zf_mi
% R1_zf_mi=RRR_zf_mi(:,:,2)./2+RRR_zf_mi(:,:,3);
% R2_zf_mi=RRR_zf_mi(:,:,2)./2+RRR_zf_mi(:,:,4);
% % RRR_zf_mp
% % RRR_mrt1
% R1_mrt1=RRR_mrt1(:,:,2)./2+RRR_mrt1(:,:,3);
% R2_mrt1=RRR_mrt1(:,:,2)./2+RRR_mrt1(:,:,4);
% % RRR_mrt2
%
%
% %%CCDF
% figure(4)
% k=6;
% %%RS
% [f_RS,x_RS]=ecdf([reshape(R1_G(:,k),[],1);reshape(R2_G(:,k),[],1)]);
% hold on;plot(x_RS,f_RS,'r--');
% %SDMA
% [f_sdma,x_sdma]=ecdf([reshape(R1_sdma(:,k),[],1);reshape(R2_sdma(:,k),[],1)]);
% plot(x_sdma,f_sdma,'b--');
% % RRR_zf
% [f_zf,x_zf]=ecdf([reshape(R1_zf(:,k),[],1);reshape(R2_zf(:,k),[],1)]);
% plot(x_zf,f_zf,'k--');
% % RRR_mrt
% [f_mrt,x_mrt]=ecdf([reshape(R1_mrt(:,k),[],1);reshape(R2_mrt(:,k),[],1)]);
% plot(x_mrt,f_mrt,'k--');
% % RRR_noma
% [f_noma,x_noma]=ecdf([reshape(R1_noma(:,k),[],1);reshape(R2_noma(:,k),[],1)]);
% plot(x_noma,f_noma,'m--');
% % RRR_zf_sp
% [f_zf_si,x_zf_si]=ecdf([reshape(R1_zf_si(:,k),[],1);reshape(R2_zf_si(:,k),[],1)]);
% plot(x_zf_si,f_zf_si,'b-.');
% % RRR_mrt1
% [f_mrt1,x_mrt1]=ecdf([reshape(R1_mrt1(:,k),[],1);reshape(R2_mrt1(:,k),[],1)]);
% plot(x_mrt1,f_mrt1,'c-.');
% % RRR_zf_mi
% [f_zf_mi,x_zf_mi]=ecdf([reshape(R1_zf_mi(:,k),[],1);reshape(R2_zf_mi(:,k),[],1)]);
% plot(x_zf_mi,f_zf_mi,'m-.');
%
% legend('RS','SDMA','RS-ZF','RS-MRT','NOMA',...
%     'SDMA-ZF','OMA-MRT1','multicasting')
% title(['SNR=',num2str(SNR(k)),', \sigma_e=',num2str(sigma_e1)]);
% xlabel(' sum rate')
% ylabel('CCDF')
%
%
% % %%%%zeroforcing
% % %%RS
% % [f_Ri,x_Ri]=ecdf(reshape(Rs__Ri(:,k),[],1));
% % [f_Rp,x_Rp]=ecdf(reshape(Rs__Rp(:,k),[],1));
% % %%NOMA
% % [f_Ni,x_Ni]=ecdf(reshape(Rs__Ni(:,k),[],1));
% % [f_Np,x_Np]=ecdf(reshape(Rs__Np(:,k),[],1));
% % %%%SDMA
% % [f_Si,x_Si]=ecdf(reshape(Rs__Si(:,k),[],1));
% % %%OMA
% % [f_Oim1,x_Oim1]=ecdf(reshape(Rs_Oim1(:,k),[],1));
% % CCDF_Oim1 = 1-f_Oim1;
% % [f_Oim,x_Oim]=ecdf(reshape(Rs_Oim(:,k),[],1));
% % %%multi
% % [f_mi,x_mi]=ecdf(reshape(Rs__mi(:,k),[],1));
% % CCDF_mi = 1-f_mi;
%
% %%%%%%%%%%%%%%%%%%%%
% % hold on
% % plot(x_RS,CCDF_RS,'r--*')
% % plot(x_Rp,CCDF_Ri,'r--o','MarkerIndices',1:round(length(CCDF_Oim)/5):length(CCDF_Oim))
% %
% %
% % hold on
% % plot(x_Ri,CCDF_Ri,'r')
% % plot(x_Rp,CCDF_Ri,'r--o','MarkerIndices',1:round(length(CCDF_Oim)/5):length(CCDF_Oim))
% % %%NOMA
% % plot(x_Ni,CCDF_Ni,'b')
% % %SDMA
% % plot(x_Si,CCDF_Si,'g')
% % %OMA
% % plot(x_Oim1,CCDF_Oim1,'c')
% % plot(x_Oim,CCDF_Oim,'c--+','MarkerIndices',1:round(length(CCDF_Oim)/5):length(CCDF_Oim))
% % %multicasting
% % plot(x_mi,CCDF_mi,'m')
% % legend('RS_{imper}','RS_{perfect}','NOMA_{imperfect}','SDMA_{imperfect}','OMA_{imperfect}','OMA_{MRTperfect}','multicasting_{imperfect}')
% % title(['\sigma_e^2=',num2str(SNR(k)),', SNR=',num2str(SNR)]);
% % xlabel(' sum rate')
% % ylabel('CCDF')



