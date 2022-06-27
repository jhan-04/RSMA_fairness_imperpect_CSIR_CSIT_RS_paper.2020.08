% clear all
% clc

Nt=2;%nuber of transmitter antenna
K=2;%number of user


SNR=[5:6:35];
%Pt=10^(SNR/10);%transmit signal power
error=sqrt(0.05);
sigma_e=error+zeros(K,1);
beta=(error)^2;
N0=1;

%gam=0;

sample=2;
noma_max=zeros(sample,6,4);
noma_max_e3=zeros(sample,6,4);
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
    
    for j=3:1:length(SNR)
        
        
        snr=SNR(j)
        Pt=10^(snr/10);
        [A,B,C,D]=cal_ABCD(Nt,K,Pt,N0,H_h,sigma_e);

        
        
        %%%%%%NOMA_two user
        [P_noma00(:,i,j),GMI_noma0(i,j,:)]=GMI_NOMA_two_user(Nt,K,H_h,Pt,N0,sigma_e);
        
        
        
        [P_noma_e,RRR_noma_e22]=GMI_NOMA_two_user(Nt,K,H_h,Pt,N0,sigma_e*0);
        [GMI_noma_e2,GMI_c_noma_e2,GMI_p_noma_e2]=cal_GMI(K,A,B,C,D,P_noma_e);
        GMI_noma_e0(i,j,:)=[GMI_noma_e2,GMI_c_noma_e2,GMI_p_noma_e2];%%%%%%%%%%%%%%%%%
        
        
        [P_noma,RRR_noma]=GMI_NOMA_two_user(Nt,K,H_h,Pt,N0,sigma_e*0);
        [GMI,GMI_c,GMI_p]=cal_GMI(K,A,B,C,D,P_noma);
        RRR_noma_e02(i,j,:)=[GMI,GMI_c,GMI_p]';%%%%%%%%%%%%%%%%%%%
        
        for k=1:100000
            amp=rand(6,1);
            pha=rand(6,1).*2.*pi;
            ran=amp.*exp(1i*pha);
            ran(3:4)=[0;0];
            
            
            pre=sqrt(Pt).*ran/norm(ran);
            
            [GMI_noma2,GMI_c_noma2,GMI_p_noma2]=cal_GMI(K,A,B,C,D, pre);
            if GMI_noma2>=noma_max(i,j,1)
                noma_max(i,j,:)= [GMI_noma2,GMI_c_noma2,GMI_p_noma2];
                prepre(:,j)=pre;
            end
            
        end
        for k=1:100000
            amp=rand(6,1);
            pha=rand(6,1).*2.*pi;
            ran=amp.*exp(1i*pha);
            ran(3:4)=[0;0];
            
            pre=sqrt(Pt).*ran/norm(ran);
            [AA,BB,CC,DD]=cal_ABCD(Nt,K,Pt,N0,H_h,sigma_e*0);
            [GMI_noma_e3,GMI_c_noma_e3,GMI_p_noma_e3]=cal_GMI(K,AA,BB,CC,DD, pre);
            
            if GMI_noma_e3>=noma_max_e3(i,j,1)
                noma_max_e3(i,j,:)= [GMI_noma_e3,GMI_c_noma_e3,GMI_p_noma_e3];
                prepre_e(:,j)=pre;
                [GMI_noma_e00,GMI_c_noma_e0,GMI_p_noma_e0]=cal_GMI(K,A,B,C,D,prepre_e(:,j));
                noma_max_e0(i,j,:)=[GMI_noma_e00,GMI_c_noma_e0,GMI_p_noma_e0];
            end
            
            
            
        end
        
    end
    
end
