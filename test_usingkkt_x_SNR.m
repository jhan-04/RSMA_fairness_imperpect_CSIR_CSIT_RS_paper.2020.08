
clear all

Nt=2;
beta=0.1;


gam_dB=-5;
gam=10^(gam_dB/20);


sample=100;

SNR=[10:6:40];




out=4;%outage rate Bound 그..머냐..동작안하는 지점...
for i=1:sample
    
    %%esimated channel
    ha=(randn(Nt,1)+1i*randn(Nt,1))/sqrt(2);%Rayleigh
    hb=(gam)*(randn(Nt,1)+1i*randn(Nt,1))/sqrt(2);
    if norm(ha)>=norm(hb)
        h1(:,i)=ha;
        h2(:,i)=hb;
    else
        h1(:,i)=hb;
        h2(:,i)=ha;
    end
    
    rho=1-abs(h1(:,i)'/norm(h1(:,i))*h2(:,i)/norm(h2(:,i)))^2;
    Rho(i)=rho;
    
    j=0;
    for snr=SNR%[0.0001,0.001,0.01,0.05,0.1,0.2]
        
        
        Pt=10^(snr/10);
        %beta=0.94*(P)^(-0.6)
        j=j+1;
        %%%%%%%%%%NRS_NOMA
        [MA_Ri(i,j),t_Ri(i,j), P1_Ri(i,j),P2_Ri(i,j), Pc_Ri(i,j),Rs_Ri(i,j),fc_Ri(:,i,j),t_Ni(i,j),Rs_Ni(i,j),fc_Ni(:,i,j)]=RSNOMA_kkt_imperfect(Pt,h1(:,i),h2(:,i),rho,beta);
        [MA_Rp(i,j),t_Rp(i,j), P1_Rp(i,j),P2_Rp(i,j), Pc_Rp(i,j),Rs_Rp(i,j),fc_Rp(:,i,j),t_Np(i,j),Rs_Np(i,j),fc_Np(:,i,j)]=RSNOMA_kkt_perfect(Pt,h1(:,i),h2(:,i),rho);
        %%%%%%
        [MA_Nim(i,j),t_Nim(i,j), P1_Nim(i,j),P2_Nim(i,j), Pc_Nim(i,j),Rs_Nim(i,j),fc_Nim(:,i,j)]=NOMA_MRT_imperfect(Pt,h1(:,i),h2(:,i),rho,beta);
        %%%%%%%%%%%%%%%%%%%%%%%%%SDMA
        [MA_Si(i,j),t_Si(i,j), P1_Si(i,j),P2_Si(i,j), Pc_Si(i,j),Rs_Si(i,j),fc_Si(:,i,j)]=SDMA_imperfect(Pt,h1(:,i),h2(:,i),rho,beta);
        [MA_Sp(i,j),t_Sp(i,j), P1_Sp(i,j),P2_Sp(i,j), Pc_Sp(i,j),Rs_Sp(i,j),fc_Sp(:,i,j)]=SDMA_perfect(Pt,h1(:,i),h2(:,i),rho);
        
        %%%%%%%OMA
        [MA_Oi(i,j),t_Oi(i,j), P1_Oi(i,j),P2_Oi(i,j), Pc_Oi(i,j),Rs_Oi(i,j),fc_Oi(:,i,j)]=OMA_imperfect(Pt,h1(:,i),h2(:,i),rho,beta);
        [MA_Op(i,j),t_Op(i,j), P1_Op(i,j),P2_Op(i,j), Pc_Op(i,j),Rs_Op(i,j),fc_Op(:,i,j)]=OMA_perfect(Pt,h1(:,i),h2(:,i),rho);
        %%%%%%%MRT
        [MA_Oim(i,j),t_Oim(i,j), P1_Oim(i,j),P2_Oim(i,j), Pc_Oim(i,j),Rs_Oim(i,j),fc_Oim(:,i,j)]=OMA_MRT_imperfect(Pt,h1(:,i),h2(:,i),rho,beta);
       
        %%%%%%multicasting
        [MA_mi(i,j),t_mi(i,j), P1_mi(i,j),P2_mi(i,j), Pc_mi(i,j),Rs_mi(i,j),fc_mi(:,i,j)]=multi_kkt_imperfect(Pt,h1(:,i),h2(:,i),rho,beta);
        [MA_mp(i,j),t_mp(i,j), P1_mp(i,j),P2_mp(i,j), Pc_mp(i,j),Rs_mp(i,j),fc_mp(:,i,j)]=multi_kkt_perfect(Pt,h1(:,i),h2(:,i),rho);
        
        
        
        %%%%%%%%%%%%%%%%%%caculate GMI
        %%%%%%%%%%RS_NOMA
        Rs__Ri(i,j)=cal_GMI(Pt,t_Ri(i,j), P1_Ri(i,j),P2_Ri(i,j), Pc_Ri(i,j),fc_Ri(:,i,j),h1(:,i),h2(:,i),rho,beta);
        Rs__Rp(i,j)=cal_GMI(Pt,t_Rp(i,j), P1_Rp(i,j),P2_Rp(i,j), Pc_Rp(i,j),fc_Rp(:,i,j),h1(:,i),h2(:,i),rho,beta);
        Rs__Ni(i,j)=cal_GMI(Pt,t_Ni(i,j), t_Ni(i,j)*Pt,0, (1-t_Ni(i,j))*Pt,fc_Ni(:,i,j),h1(:,i),h2(:,i),rho,beta);
        Rs__Np(i,j)=cal_GMI(Pt,t_Np(i,j), t_Np(i,j)*Pt,0, (1-t_Np(i,j))*Pt,fc_Np(:,i,j),h1(:,i),h2(:,i),rho,beta);
        %MRT
        %%Rs__Nim(i,j)=cal_GMI(P,t_Nim(i,j), t_Nim(i,j)*P,0, (1-t_Nim(i,j))*P,fc_Nim(:,i,j),h1(:,i),h2(:,i),rho,beta);
        %%%%%%%%%%%%%%%%%%%%%%%%%SDMA
        Rs__Si(i,j)=cal_GMI(Pt,t_Si(i,j), P1_Si(i,j),P2_Si(i,j), Pc_Si(i,j),fc_Si(:,i,j),h1(:,i),h2(:,i),rho,beta);
        Rs__Sp(i,j)=cal_GMI(Pt,t_Sp(i,j), P1_Sp(i,j),P2_Sp(i,j), Pc_Sp(i,j),fc_Sp(:,i,j),h1(:,i),h2(:,i),rho,beta);
        
        %%%%%%%%%%%%%%%%%%%OMA
        Rs__Oi(i,j)=cal_GMI(Pt,t_Oi(i,j), P1_Oi(i,j),P2_Oi(i,j), Pc_Oi(i,j),fc_Oi(:,i,j),h1(:,i),h2(:,i),rho,beta);
        Rs__Op(i,j)=cal_GMI(Pt,t_Op(i,j), P1_Op(i,j),P2_Op(i,j), Pc_Op(i,j),fc_Op(:,i,j),h1(:,i),h2(:,i),rho,beta);
        
        %%%%%%%%%%%%%%%multi
        Rs__mi(i,j)=cal_GMI(Pt,t_mi(i,j), P1_mi(i,j),P2_mi(i,j), Pc_mi(i,j),fc_mi(:,i,j),h1(:,i),h2(:,i),rho,beta);
        Rs__mp(i,j)=cal_GMI(Pt,t_mp(i,j), P1_mp(i,j),P2_mp(i,j), Pc_mp(i,j),fc_mp(:,i,j),h1(:,i),h2(:,i),rho,beta);
        
        
        
                fprintf('=======sample=%d===================SNR=%d===================\n\n',i,snr)
        %         fprintf('RS_KKT_im: MA=%d , t=%d, P1=%d, P2=%d Pc=%d Rs=%d GMI=%d \n ',MA_Ri(i,j),t_Ri(i,j), P1_Ri(i,j),P2_Ri(i,j), Pc_Ri(i,j),Rs_Ri(i,j),Rs__Ri(i,j))
        %         fprintf('NOMA_KKT_im: t=%d, Rs=%d Rs2=%d\n\n',t_Ni(i,j),Rs_Ni(i,j),Rs__Ni(i,j))
        %         fprintf('RS_KKT: MA=%d , t=%d, P1=%d, P2=%d Pc=%d Rs=%d GMI=%d\n ',MA_Rp(i,j),t_Rp(i,j), P1_Rp(i,j),P2_Rp(i,j), Pc_Rp(i,j),Rs_Rp(i,j),Rs__Rp(i,j))
        %         fprintf('NOMA_KKT: t=%d, Rs=%d GMI=%d\n\n',t_Np(i,j),Rs_Np(i,j),Rs__Np(i,j))
        %
        %         fprintf('SDMA_im: MA=%d , t=%d, P1=%d, P2=%d Pc=%d Rs=%d GMI=%d\n ',MA_Si(i,j),t_Si(i,j), P1_Si(i,j),P2_Si(i,j), Pc_Si(i,j),Rs_Si(i,j),Rs__Si(i,j))
        %         fprintf('SDMA: MA=%d , t=%d, P1=%d, P2=%d Pc=%d Rs=%d GMI=%d\n\n ',MA_Sp(i,j),t_Sp(i,j), P1_Sp(i,j),P2_Sp(i,j), Pc_Sp(i,j),Rs_Sp(i,j),Rs__Sp(i,j))
        %
        %         fprintf('OMA_im: MA=%d , t=%d, P1=%d, P2=%d Pc=%d Rs=%d GMI=%d \n ',MA_Oi(i,j),t_Oi(i,j), P1_Oi(i,j),P2_Oi(i,j), Pc_Oi(i,j),Rs_Oi(i,j),Rs__Oi(i,j))
        %         fprintf('OMA: MA=%d , t=%d, P1=%d, P2=%d Pc=%d Rs=%d GMI=%d\n\n ',MA_Op(i,j),t_Op(i,j), P1_Op(i,j),P2_Op(i,j), Pc_Op(i,j),Rs_Op(i,j),Rs__Op(i,j))
        %
        %         fprintf('multi_KKT_im: MA=%d , t=%d, P1=%d, P2=%d Pc=%d Rs=%d GMI=%d\n ',MA_mi(i,j),t_mi(i,j), P1_mi(i,j),P2_mi(i,j), Pc_mi(i,j),Rs_mi(i,j),Rs__mi(i,j))
        %         fprintf('multi_KKT: MA=%d , t=%d, P1=%d, P2=%d Pc=%d Rs=%d GMI=%d\n\n ',MA_mp(i,j),t_mp(i,j), P1_mp(i,j),P2_mp(i,j), Pc_mp(i,j),Rs_mp(i,j),Rs__mp(i,j))
        
        %in inperfect optimization progress, Rs shoud be same with GMI(=>Rs_Ri==Rs__Ri)
        
    end
    
end

%%%number of MA cases in samples
figure(1)
for k=1:1:numel(SNR)
    
    SDMA_case(k)=numel(find (MA_Ri(:,k)==1));
    NOMA_case(k)=numel(find (MA_Ri(:,k)==2));
    OMA_case(k)=numel(find (MA_Ri(:,k)==3));
    mul_case(k)=numel(find (MA_Ri(:,k)==4));
    RS_case(k)=numel(find (MA_Ri(:,k)==5));
    %SDMA_case(k)+NOMA_case(k)+OMA_case(k)+mul_case(k)+RS_case(k)==sample;
end
CASE=[SDMA_case;NOMA_case;OMA_case;mul_case;RS_case]';
bar(SNR,CASE,'stacked')
legend('SDMA','NOMA','OMA','multicasting','RS')
xlabel('SNR')
title(['beta(channel error covariance parameter)=',num2str(beta),', channel disparty=',num2str(gam_dB)])

%%%%average of GMI with sample channel
figure(2)
hold on
%%RS
plot(SNR,mean(Rs__Ri),'r')
plot(SNR,mean(Rs__Rp),'r--o')
%%NOMA
plot(SNR,mean(Rs__Ni),'b')
plot(SNR,mean(Rs__Np),'b--o')
plot(SNR,mean(Rs_Nim),'b--d')
%SDMA
plot(SNR,mean(Rs__Si),'g')
plot(SNR,mean(Rs__Sp),'g--o')
%OMA
plot(SNR,mean(Rs__Oi),'c')
plot(SNR,mean(Rs__Op),'c--o')
%MRT
plot(SNR,mean(Rs_Oim),'c--+')
%multicasting
plot(SNR,mean(Rs__mi),'m')
plot(SNR,mean(Rs__mp),'m--o')
legend('RS_{imper}','RS_{perfect}','NOMA_{imperfect}','NOMA_{perfect}','NOMA_{MRTimperfect}','SDMA_{imperfect}','SDMA_{perfect}','OMA_{imperfect}','OMA_{perfect}','OMA_{MRTperfect}','multicasting_{imperfect}','multicasting_{perfect}')
title(['beta(channel error covariance parameter)=',num2str(beta),', channel disparty=',num2str(gam_dB)])
xlabel('SNR')
ylabel('GMI')
hold off





figure(3)

%outage
OUT_Ri=sum(Rs__Ri>out);
OUT_Rp=sum(Rs__Rp>out);

OUT_Ni=sum(Rs__Ni>out);
OUT_Np=sum(Rs__Np>out);
OUT_Nim=sum(Rs_Nim>out);

OUT_Si=sum(Rs__Si>out);
OUT_Sp=sum(Rs__Sp>out);

OUT_Oi=sum(Rs__Oi>out);
OUT_Op=sum(Rs__Op>out);
OUT_Oim=sum(Rs_Oim>out);

OUT_mi=sum(Rs__mi>out);
OUT_mp=sum(Rs__mp>out);
OUT=[OUT_Ri;OUT_Rp;OUT_Ni;OUT_Np;OUT_Nim;OUT_Si;OUT_Sp;OUT_Oi;OUT_Op;OUT_Oim;OUT_mi;OUT_mp]';
bar(SNR,OUT)

legend('RS_{imper}','RS_{perfect}','NOMA_{imperfect}','NOMA_{perfect}','NOMA_{MRTimperfect}','SDMA_{imperfect}','SDMA_{perfect}','OMA_{imperfect}','OMA_{perfect}','OMA_{MRTperfect}','multicasting_{imperfect}','multicasting_{perfect}')
title(['beta(channel error covariance parameter)=',num2str(beta),', channel disparty=',num2str(gam_dB)])
xlabel('SNR')
ylabel('GMI')


%%CCDF
figure(4)
k=2;
%%RS
[f_Ri,x_Ri]=ecdf(reshape(Rs__Ri(:,k),[],1));
CCDF_Ri = 1-f_Ri;
[f_Rp,x_Rp]=ecdf(reshape(Rs__Rp(:,k),[],1));
CCDF_Ri = 1-f_Rp;
%%NOMA
[f_Ni,x_Ni]=ecdf(reshape(Rs__Ni(:,k),[],1));
CCDF_Ni = 1-f_Ni;
[f_Nim,x_Nim]=ecdf(reshape(Rs_Nim(:,k),[],1));
CCDF_Nim = 1-f_Nim;
%%%SDMA
[f_Si,x_Si]=ecdf(reshape(Rs__Si(:,k),[],1));
CCDF_Si = 1-f_Si;
%%OMA
[f_Oi,x_Oi]=ecdf(reshape(Rs__Oi(:,k),[],1));
CCDF_Oi = 1-f_Oi;
[f_Oim,x_Oim]=ecdf(reshape(Rs_Oim(:,k),[],1));
CCDF_Oim = 1-f_Oim;
%%multi
[f_mi,x_mi]=ecdf(reshape(Rs__mi(:,k),[],1));
CCDF_mi = 1-f_mi;

%%%%%%%%%%%%%%%%%%%%
hold on
plot(x_Ri,CCDF_Ri,'r')
plot(x_Rp,CCDF_Ri,'r--o','MarkerIndices',1:round(length(CCDF_Oim)/5):length(CCDF_Oim))
%%NOMA
plot(x_Ni,CCDF_Ni,'b')
plot(x_Nim,CCDF_Nim,'b--d','MarkerIndices',1:round(length(CCDF_Oim)/5):length(CCDF_Oim))
%SDMA
plot(x_Si,CCDF_Si,'g')
%OMA
plot(x_Oi,CCDF_Oi,'c')
plot(x_Oim,CCDF_Oim,'c--+','MarkerIndices',1:round(length(CCDF_Oim)/5):length(CCDF_Oim))
%multicasting
plot(x_mi,CCDF_mi,'m')
legend('RS_{imper}','RS_{perfect}','NOMA_{imperfect}','NOMA_{MRTimperfect}','SDMA_{imperfect}','OMA_{imperfect}','OMA_{MRTperfect}','multicasting_{imperfect}')
title(['beta(channel error covariance parameter)=',num2str(beta),', channel disparty=',num2str(gam_dB),', SNR=',num2str(SNR(k))]);
xlabel(' sum rate')
ylabel('CCDF')
%     Rs__Ri
%     Rs__Ni
%     Rs__Rp
%     Rs__Np
%     Rs__Si
%     Rs__Sp
%     Rs__Oi
%     Rs__Op
%     Rs__mi
%     Rs__mp












