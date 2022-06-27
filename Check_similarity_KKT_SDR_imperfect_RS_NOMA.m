clc


Nt=2;

P=100;

gam_dB=(-rand()+0.5)*10
gam=10^(gam_dB/20);
beta=0.1;

h1=(randn(Nt,1)+1i*randn(Nt,1))/sqrt(2);%Rayleigh
h2=(gam)*(randn(Nt,1)+1i*randn(Nt,1))/sqrt(2);


rho=1-abs(h1'/norm(h1)*h2/norm(h2))^2;
%%%%%%%%%%%%%%%%%%%%%%%

[MA_x,t_x, P1_x,P2_x, Pc_x,Rs_x,fc_x,t2,Rs_2,fc_2]=RSNOMA_sdr_imperfect(P,h1,h2,rho,beta);
fprintf('RS_SDR: MA=%d , t=%d, P1=%d, P2=%d Pc=%d Rs=%d \n ',MA_x,t_x, P1_x,P2_x, Pc_x,Rs_x) 
fprintf('NOMA_SDR: t=%d, Rs=%d \n ',t2,Rs_2)

[MA_x,t_x, P1_x,P2_x, Pc_x,Rs_x,fc_x]=NOMA_sdr_imperfect(P,h1,h2,rho,beta);
fprintf('NOMA_SDR2: t=%d, Rs=%d \n ',t_x,Rs_x)


%%%%%%%%%%%%%%
[MA_x,t_x, P1_x,P2_x, Pc_x,Rs_x,fc_x,t2,Rs_2,fc_2]=RSNOMA_kkt_imperfect(P,h1,h2,rho,beta);
fprintf('RS_KKT: MA=%d , t=%d, P1=%d, P2=%d Pc=%d Rs=%d \n ',MA_x,t_x, P1_x,P2_x, Pc_x,Rs_x) 
fprintf('NOMA_KKT: t=%d, Rs=%d \n ',t2,Rs_2)


[MA_x,t_x, P1_x,P2_x, Pc_x,Rs_x,fc_x]=NOMA_kkt_imperfect(P,h1,h2,rho,beta);
fprintf('NOMA_KKT2: t=%d, Rs=%d \n ',t_x,Rs_x)












