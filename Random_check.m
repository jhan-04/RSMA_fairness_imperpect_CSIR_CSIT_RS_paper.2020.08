
GMI_ran=zeros(sample,length(SNR),K+2);
p_ran_final=zeros(Nt*(K+1),1);
C_ran_final=zeros(K,1);

for r=1:10^6
    
    if mod(r,2000)==0
        fprintf('random= %d  \n',r)
    end
    for i=1:sample
        
        i;
        %user1 strong channel
        %user2 week channel
        
        
        H_h(:,1)=h1(:,i);
        H_h(:,2)=h2(:,i);
        
        for j=1:1:length(SNR)
            
            
            snr=SNR(j);
            Pt=10^(snr/10);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Optimization GMI
            
            
            p_amp=rand(Nt*(K+1),1);
            p_pha=2.*pi.*rand(Nt*(K+1),1);
            p_ran0= p_amp.*exp(1i*p_pha);
            
            p_ran=sqrt(Pt)*p_ran0/norm(p_ran0);
            
            C_rand=rand(K,1);
            
            [GMI_min,GMI_c,GMI_p]=cal_GMI_min(K,A,B,C,D,p_ran,CC);
            
            if GMI_min>GMI_ran(i,j,:)
                GMI_ran(i,j,:)=[GMI_min,GMI_c,GMI_p]';
                p_ran_final(:,i,j)=p_ran;
                C_ran_final(:,i,j)=C_rand;
            end
            
            
        end
        

        
        
        
        
    end
    %norm(h1(:,i))/norm(h2(:,i))
end

%figure
hold on

plot(SNR,mean(RRR_RS(:,:,1)),':hk','LineWidth',1,'MarkerSize',3);


ylim[1,3.5]