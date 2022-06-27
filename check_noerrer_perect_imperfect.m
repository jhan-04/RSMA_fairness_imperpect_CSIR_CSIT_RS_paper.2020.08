


Nt=2;
i=1;
j=1;
P=1000;

gam_dB=-5
gam=10^(gam_dB/20);
beta=0.1;

% h1=(randn(Nt,1)+1i*randn(Nt,1))/sqrt(2);%Rayleigh
% h2=(gam)*(randn(Nt,1)+1i*randn(Nt,1))/sqrt(2);


rho=1-abs(h1(:,i)'/norm(h1(:,i))*h2(:,i)/norm(h2(:,i)))^2;
%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
 [MA_Ri(i,j),t_Ri(i,j), P1_Ri(i,j),P2_Ri(i,j), Pc_Ri(i,j),Rs_Ri(i,j),fc_Ri(:,i,j),t_Ni(i,j),Rs_Ni(i,j),fc_Ni(:,i,j)]=RSNOMA_sdr_imperfect(P,h1(:,i),h2(:,i),rho,beta);       
%[MA_x,t_x, P1_x,P2_x, Pc_x,Rs_x,fc_x,t2,Rs_2,fc_2]=RSNOMA_sdr_imperfect(P,h1,h2,rho,beta);
fprintf('RS_SDR_im: MA=%d , t=%d, P1=%d, P2=%d Pc=%d Rs=%d \n ',MA_Ri(i,j),t_Ri(i,j), P1_Ri(i,j),P2_Ri(i,j), Pc_Ri(i,j),Rs_Ri(i,j)) 
fprintf('NOMA_SDR_im: t=%d, Rs=%d \n ',t_Ni(i,j),Rs_Ni(i,j))

%%%%%%%%%%%%%%%%%% 

        [MA_Rp(i,j),t_Rp(i,j), P1_Rp(i,j),P2_Rp(i,j), Pc_Rp(i,j),Rs_Rp(i,j),fc_Rp(:,i,j),t_Np(i,j),Rs_Np(i,j),fc_Np(:,i,j)]=RSNOMA_sdr_perfect(P,h1(:,i),h2(:,i),rho);
fprintf('RS_SDR: MA=%d , t=%d, P1=%d, P2=%d Pc=%d Rs=%d \n ',MA_Rp(i,j),t_Rp(i,j), P1_Rp(i,j),P2_Rp(i,j), Pc_Rp(i,j),Rs_Rp(i,j)) 
fprintf('NOMA_SDR_: t=%d, Rs=%d \n ',t_Np(i,j),Rs_Np(i,j))

% [MA_x,t_x, P1_x,P2_x, Pc_x,Rs_x,fc_x]=RS_sdr_perfect(P,h1,h2,rho);
% fprintf('RS_SDR2: t=%d, Rs=%d \n ',t_x,Rs_x)
% [MA_x,t_x, P1_x,P2_x, Pc_x,Rs_x,fc_x]=NOMA_sdr_perfect(P,h1,h2,rho);
% fprintf('NOMA_SDR2: t=%d, Rs=%d \n ',t_x,Rs_x)

fprintf('\n===============================================\n')

%%%%%%%%%%%%%%
[MA_Rik(i,j),t_Rik(i,j), P1_Rik(i,j),P2_Rik(i,j), Pc_Rik(i,j),Rs_Rik(i,j),fc_Rik(:,i,j),t_Nik(i,j),Rs_Nik(i,j),fc_Nik(:,i,j)]=RSNOMA_kkt_imperfect(P,h1(:,i),h2(:,i),rho,beta);
%[MA_x,t_x, P1_x,P2_x, Pc_x,Rs_x,fc_x,t2,Rs_2,fc_2]=RSNOMA_kkt_imperfect(P,h1,h2,rho,beta);
fprintf('RS_KKT_im: MA=%d , t=%d, P1=%d, P2=%d Pc=%d Rs=%d \n ',MA_Rik(i,j),t_Rik(i,j), P1_Rik(i,j),P2_Rik(i,j), Pc_Rik(i,j),Rs_Rik(i,j)) 
fprintf('NOMA_KKT_im: t=%d, Rs=%d \n ',t_Nik(i,j),Rs_Nik(i,j))


%%%%%%%%%%%%%%%%%%%


        [MA_Rpk(i,j),t_Rpk(i,j), P1_Rpk(i,j),P2_Rpk(i,j), Pc_Rpk(i,j),Rs_Rpk(i,j),fc_Rpk(:,i,j),t_Npk(i,j),Rs_Npk(i,j),fc_Npk(:,i,j)]=RSNOMA_kkt_perfect(P,h1(:,i),h2(:,i),rho);
        
%[MA_x,t_x, P1_x,P2_x, Pc_x,Rs_x,fc_x,t2,Rs_2,fc_2]=RSNOMA_kkt_perfect(P,h1,h2,rho);
fprintf('RS_kkt: MA=%d , t=%d, P1=%d, P2=%d Pc=%d Rs=%d \n ',MA_Rpk(i,j),t_Rpk(i,j), P1_Rpk(i,j),P2_Rpk(i,j), Pc_Rpk(i,j),Rs_Rpk(i,j)) 
fprintf('NOMA_kkt: t=%d, Rs=%d \n ',t_Npk(i,j),Rs_Npk(i,j))

% [MA_x,t_x, P1_x,P2_x, Pc_x,Rs_x,fc_x]=RS_kkt_perfect(P,h1,h2,rho);
% 
% fprintf('RS_kkt2: t=%d, Rs=%d \n ',t_x,Rs_x)
% [MA_x,t_x, P1_x,P2_x, Pc_x,Rs_x,fc_x]=NOMA_kkt_perfect(P,h1,h2,rho);
% fprintf('NOMA_kkt2: t=%d, Rs=%d \n ',t_x,Rs_x)














