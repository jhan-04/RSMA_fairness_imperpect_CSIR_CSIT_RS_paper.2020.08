clear all
clc


Nt=5;%nuber of transmitter antenna
M=5;%number of user
SNR=10;
Pt=10^(SNR/10);%transmit signal power
N0=1;%gaussina noise

sigma_e1=[0.1,0.2:0.2:0.6]; %sigma_e2=0.3; sigma_e3=0.3;
sigma_e=randn(M,1)*0+sigma_e1;%M X 1 vector
result_set_e=zeros(40);
result_nc_set_e=zeros(40);
result_zf_set_e=zeros(40);
sample=50;
%for j=1:1:length(sigma_e1)

for i=1:sample
    
    i
    %%%hK
    for k=1:M
        h1(:,k)=(randn(Nt,1)+1i*randn(Nt,1))/sqrt(2);%Rayleigh h1
    end
    
    
    for j=1:1:length(sigma_e1)
        sigma_e=randn(M,1)*0+sigma_e1(j);%M X 1 vector
        for k=1:1:M
            h(:,k)=sqrt(1-sigma_e(k)^2)*h1(:,k);
        end
        
        [A,B,C,D]=cal_ABCD(Nt,M,Pt,N0,h,sigma_e);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Optimization GMI
        
        %%%%%%%%%%%%Rate splitting%%%%%%%%%%%%%%
        [X,result]=cal_X(Nt,M,Pt,A,B,C,D);
        result_set(1:length(result),i,j)=result;
        %reult/log(2)==GMI
        [GMI_X(i,j)]=cal_GMI_withX(M,A,B,C,D,X);
        [p_max(:,i,j)]=find_p(M,Pt,A,B,C,D,X);norm(p_max(:,i,j))^2;
        [GMI_L(i,j)]=cal_GMI(M,A,B,C,D,p_max(:,i,j));
        %         result(length(result))/log(2)
        %         GMI_L(i,j)
        
        
        %%%%%%%%%%%%no_common_part(=SDMA)%%%%%%%%%%%%%%
        [X_nc,result_nc]=cal_X_no_common(Nt,M,Pt,A,B,C,D);
        result_nc_set(1:length(result_nc),i,j)=result_nc;
        [GMI_X_nc(i,j)]=cal_GMI_withX(M,A,B,C,D,X_nc);
        [p_nc(:,i,j)]=find_p(M,Pt,A,B,C,D,X_nc);norm(p_nc(:,i,j))^2;
        [GMI_L_nc(i,j)]=cal_GMI(M,A,B,C,D,p_nc(:,i,j));
        %                 result_nc(length(result_nc))/log(2)
        %                 GMI_L_nc(i,j)
        
        
        %%%%%%%%%Zero_Forcing%%%%%%%%%%
        [Xc,PK,result_zf]=cal_X_ZF(Nt,M,Pt,h,N0,sigma_e);
        result_zf_set(1:length(result_zf),i,j)=result_zf;
        pk_zf=h*inv(h'*h);
        for k=1:M
            mag(k)=norm(pk_zf(:,k));
        end
        pk(:,i,j)=reshape(pk_zf.*(sqrt(PK)./mag')',[],1);
        [p_zf(:,i,j)]=find_p_ZF(M,Pt,A,B,C,D,Xc,pk(:,i,j));norm(p_zf(:,i,j))^2;
        [GMI_zf(i,j)]=cal_GMI(M,A,B,C,D,p_zf(:,i,j));
        %                 result_zf(length(result_zf))/log(2)
        %                 GMI_zf(i,j)
        
        
        
        
        %%%%%%%%%MRT%%%%%%%%%%
        [Xc,PK,result_mrt]=cal_X_MRT(Nt,M,Pt,h,N0,sigma_e);
        result_mrt_set(1:length(result_mrt),i,j)=result_mrt;
        for k=1:M
            mag(k)=norm(h(:,k));
        end
        pk(:,i,j)=reshape(h.*(sqrt(PK)./mag')',[],1);
        [p_mrt(:,i,j)]=find_p_ZF(M,Pt,A,B,C,D,Xc,pk(:,i,j));norm(p_mrt(:,i,j))^2;
        [GMI_mrt(i,j)]=cal_GMI(M,A,B,C,D,p_mrt(:,i,j));
        %                 result_mrt(length(result_mrt))/log(2)
        %                 GMI_mrt(i,j)
        
        
        
        
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%Caculation achevable Rate
        for m=1:500
            for k=1:M
                hh(:,k)=h(:,k)+sigma_e(k)*(randn(Nt,1)+1i*randn(Nt,1))/sqrt(2);
            end
            %%%%%%%%%%%%%%Rate splitting%%%%%%%%%%%%%%%%
            [Rs_set(:,i,j)]=cal_ach_rate(Nt,M,hh,N0,p_max(:,i,j));
            ach_rate_set(i,j,m)=sum(Rs_set(:,i,j));
            
            %%%%%%%%%%%%%%no_common_part(=SDMA)%%%%%%%%%%%%%%
            [Rs_set_nc(:,i,j)]=cal_ach_rate(Nt,M,hh,N0,p_nc(:,i,j));
            ach_rate_nc_set(i,j,m)=sum(Rs_set_nc(:,i,j));
            
            %%%%%%%%%Zero_Forcing%%%%%%%%%%
            [Rs_set_zf(:,i,j)]=cal_ach_rate(Nt,M,hh,N0,p_zf(:,i,j));
            ach_rate_zf_set(i,j,m)=sum(Rs_set_zf(:,i,j));
            
            %%%%%%%%%MRT%%%%%%%%%%
            [Rs_set_mrt(:,i,j)]=cal_ach_rate(Nt,M,hh,N0,p_mrt(:,i,j));
            ach_rate_mrt_set(i,j,m)=sum(Rs_set_mrt(:,i,j));
            
        end
        
        
        ach_rate(i,j)=mean(ach_rate_set(i,j,:));
        ach_rate_nc(i,j)=mean(ach_rate_nc_set(i,j,:));
        ach_rate_zf(i,j)=mean(ach_rate_zf_set(i,j,:));
        ach_rate_mrt(i,j)=mean(ach_rate_mrt_set(i,j,:));
        
        fprintf('================sample=%d=======sigma=%d=============\n',i,sigma_e1(j))
        fprintf('RS: GMI=%d, ach_rate=%d,\n NO_RS: GMI=%d, ach_rate=%d \n',GMI_L(i,j),ach_rate(i,j),GMI_X_nc(i,j),ach_rate_nc(i,j))
        fprintf(' RS_ZF: GMI=%d, ach_rate=%d\n  RS_MRT: GMI=%d, ach_rate=%d\n',GMI_zf(i,j),ach_rate_zf(i,j),GMI_mrt(i,j),ach_rate_mrt(i,j))
        
        %
        %         fprintf('RS: GMI=%d, \n NO_RS: GMI=%d \n',GMI_L(i,j),GMI_X_nc(i,j))
        %         fprintf(' RS_ZF: GMI=%d\n  RS_MRT: GMI=%d\n',GMI_zf(i,j),GMI_mrt(i,j))
        
        
        
        %===========no information about channel error============================
        [AA,BB,CC,DD]=cal_ABCD(Nt,M,Pt,N0,h,sigma_e*0);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Optimization GMI
        
        %%%%%%%%%%%%Rate splitting%%%%%%%%%%%%%%
        [X_e,result_e]=cal_X(Nt,M,Pt,AA,BB,CC,DD);
        result_set_e(1:length(result_e),i,j)=result_e;
        %reult/log(2)==GMI
        [GMI_X_e(i,j)]=cal_GMI_withX(M,AA,BB,CC,DD,X_e);
        [p_max_e(:,i,j)]=find_p(M,Pt,AA,BB,CC,DD,X_e);
        [GMI_L_e(i,j)]=cal_GMI(M,AA,BB,CC,DD,p_max_e(:,i,j));
        %         result(length(result))/log(2)
        %         GMI_L(i,j)
        
        
        %%%%%%%%%%%%no_common_part(=SDMA)%%%%%%%%%%%%%%
        [X_nc_e,result_nc_e]=cal_X_no_common(Nt,M,Pt,AA,BB,CC,DD);
        result_nc_set_e(1:length(result_nc_e),i,j)=result_nc_e;
        [GMI_X_nc_e(i,j)]=cal_GMI_withX(M,AA,BB,CC,DD,X_nc_e);
        [p_nc_e(:,i,j)]=find_p(M,Pt,AA,BB,CC,DD,X_nc_e);
        [GMI_L_nc_e(i,j)]=cal_GMI(M,AA,BB,CC,DD,p_nc_e(:,i,j));
        %         result_nc(length(result_nc))/log(2)
        %         GMI_L_nc(i,j)
        
        
        %%%%%%%%%Zero_Forcing%%%%%%%%%%
        [Xc_e,PK_e,result_zf_e]=cal_X_ZF(Nt,M,Pt,h,N0,sigma_e*0);
        result_zf_set_e(1:length(result_zf_e),i,j)=result_zf_e;
        pk_zf_e=h*inv(h'*h);
        for k=1:M
            mag(k)=norm(pk_zf_e(:,k));
        end
        pk_e(:,i,j)=reshape(pk_zf_e.*(sqrt(PK_e)./mag')',[],1);
        [p_zf_e(:,i,j)]=find_p_ZF(M,Pt,AA,BB,CC,DD,Xc_e,pk_e(:,i,j));
        [GMI_zf_e(i,j)]=cal_GMI(M,AA,BB,CC,DD,p_zf_e(:,i,j));
        %         result_zf(length(result_zf))/log(2)
        %         GMI_zf(i,j)
        
        
        %%%%%%%%%MRT%%%%%%%%%%
        [Xc_e,PK_e,result_mrt_e]=cal_X_MRT(Nt,M,Pt,h,N0,sigma_e*0);
        result_mrt_set_e(1:length(result_mrt_e),i,j)=result_mrt_e;
        for k=1:M
            mag(k)=norm(h(:,k));
        end
        pk(:,i,j)=reshape(h.*(sqrt(PK)./mag')',[],1);
        [p_mrt_e(:,i,j)]=find_p_ZF(M,Pt,A,B,C,D,Xc,pk(:,i,j));
        [GMI_mrt_e(i,j)]=cal_GMI(M,A,B,C,D,p_mrt_e(:,i,j));
        %                 result_mrt_e(length(result_mrt_e))/log(2)
        %                 GMI_mrt_e(i,j)
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%Caculation achevable Rate
        for m=1:500
            for k=1:M
                hh(:,k)=h(:,k)+sigma_e(k)*(randn(Nt,1)+1i*randn(Nt,1))/sqrt(2);
            end
            %%%%%%%%%%%%%%Rate splitting%%%%%%%%%%%%%%%%
            [Rs_set_e(:,i,j)]=cal_ach_rate(Nt,M,hh,N0,p_max_e(:,i,j));
            ach_rate_set_e(i,j,m)=sum(Rs_set_e(:,i,j));
            
            %%%%%%%%%%%%%%no_common_part(=SDMA)%%%%%%%%%%%%%%
            [Rs_set_nc_e(:,i,j)]=cal_ach_rate(Nt,M,hh,N0,p_nc_e(:,i,j));
            ach_rate_nc_set_e(i,j,m)=sum(Rs_set_nc_e(:,i,j));
            
            %%%%%%%%%Zero_Forcing%%%%%%%%%%
            [Rs_set_zf_e(:,i,j)]=cal_ach_rate(Nt,M,hh,N0,p_zf_e(:,i,j));
            ach_rate_zf_set_e(i,j,m)=sum(Rs_set_zf_e(:,i,j));
            %%%%%%%%%MRT%%%%%%%%%%
            [Rs_set_mrt_e(:,i,j)]=cal_ach_rate(Nt,M,hh,N0,p_mrt_e(:,i,j));
            ach_rate_mrt_set_e(i,j,m)=sum(Rs_set_mrt_e(:,i,j));
            
        end
        
        
        ach_rate_e(i,j)=mean(ach_rate_set_e(i,j,:));
        ach_rate_nc_e(i,j)=mean(ach_rate_nc_set_e(i,j,:));
        ach_rate_zf_e(i,j)=mean(ach_rate_zf_set_e(i,j,:));
        ach_rate_mrt_e(i,j)=mean(ach_rate_mrt_set_e(i,j,:));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf('=======no error information======sample=%d=======sigma=%d=============\n',i,sigma_e1(j))
        fprintf('RS: GMI=%d, ach_rate=%d,\n NO_RS: GMI=%d, ach_rate=%d \n',GMI_L_e(i,j),ach_rate_e(i,j),GMI_X_nc_e(i,j),ach_rate_nc_e(i,j))
        fprintf(' RS_ZF: GMI=%d, ach_rate=%d\n  RS_MRT: GMI=%d, ach_rate=%d\n',GMI_zf_e(i,j),ach_rate_zf_e(i,j),GMI_mrt_e(i,j),ach_rate_mrt_e(i,j))
        
        %
        
        
        
        
        
    end
    
    
    
    
    
end
%i=1:sample;

%i=1:sample;
hold on
plot(sigma_e1,mean(GMI_L),'-.r')
plot(sigma_e1,mean(GMI_L_nc),'b-.')
plot(sigma_e1,mean(GMI_zf),'g-.')
plot(sigma_e1,mean(GMI_mrt),'c-.')

plot(sigma_e1,mean(ach_rate),'r-o')
plot(sigma_e1,mean(ach_rate_nc),'b-o')
plot(sigma_e1,mean(ach_rate_zf),'g-o')
plot(sigma_e1,mean(ach_rate_zf),'c-o')

%plot(SNR,mean(GMI_L_e),'-r',SNR,mean(GMI_L_nc_e),'b-',SNR,mean(GMI_zf_e),'g-')
plot(sigma_e1,mean(ach_rate_e),'--r>')
plot(sigma_e1,mean(ach_rate_nc_e),'b--^')
plot(sigma_e1,mean(ach_rate_zf_e),'g--<')
plot(sigma_e1,mean(ach_rate_mrt_e),'c--<')
% fprintf('=======sample=====================================\n\n')
% fprintf('RS_KKT_im: MA=%d , t=%d, P1=%d, P2=%d Pc=%d Rs=%d GMI=%d \n ',)

legend('GMI RS','GMI no-RS','GMI RS-ZF','GMI RS-MRT','rate RS','rate no-RS','rate RS-ZF','rate RS-MRT','rate RS no info','rate no-RS  no info','rate RS-ZF no info','rate RS-MRT no info')
xlabel('sigma_e')
ylabel('rate (bits/s/Hz)')
title(['Nt=',num2str(Nt),', M=',num2str(M),', SNR=',num2str(SNR),', sample=',num2str(sample)])

figure(2)
hold on
plot(sigma_e1,mean(GMI_L_e),'-.r')
plot(sigma_e1,mean(GMI_L_nc_e),'b-.')
plot(sigma_e1,mean(GMI_zf_e),'g-.')
plot(sigma_e1,mean(GMI_mrt_e),'c-.')
legend('GMI RS','GMI no-RS','GMI RS-ZF','GMI RS-MRT')
xlabel('sigma_e')
ylabel('rate (bits/s/Hz)')
title(['Nt=',num2str(Nt),', M=',num2str(M),', SNR=',num2str(SNR),', sample=',num2str(sample)])

