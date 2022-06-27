
Nt=2;

P=1000;

gam_dB=-(rand()+0.5)*10;
gam=10^(gam_dB/20);
beta=0;

h1=(randn(Nt,1)+1i*randn(Nt,1))/sqrt(2);%Rayleigh
h2=(gam)*(randn(Nt,1)+1i*randn(Nt,1))/sqrt(2);


rho=1-abs(h1'/norm(h1)*h2/norm(h2))^2


fprintf('\n==============START=============\n\n')

[MA_x,t_x, P1_x,P2_x, Pc_x,Rs_x,fc_x]=SDMA_perfect(P,h1,h2,rho);
fprintf('SDMA: MA=%d , t=%d, P1=%d, P2=%d Pc=%d Rs=%d \n ',MA_x,t_x, P1_x,P2_x, Pc_x,Rs_x) 
[MA_x,t_x, P1_x,P2_x, Pc_x,Rs_x,fc_x]=SDMA_imperfect(P,h1,h2,rho,beta);
fprintf('SDMA_im: MA=%d , t=%d, P1=%d, P2=%d Pc=%d Rs=%d \n ',MA_x,t_x, P1_x,P2_x, Pc_x,Rs_x) 



fprintf('\n===========================\n\n')

[MA_x,t_x, P1_x,P2_x, Pc_x,Rs_x,fc_x]=multi_sdr_perfect(P,h1,h2,rho);
fprintf('multi: MA=%d , t=%d, P1=%d, P2=%d Pc=%d Rs=%d \n ',MA_x,t_x, P1_x,P2_x, Pc_x,Rs_x)
[MA_x,t_x, P1_x,P2_x, Pc_x,Rs_x,fc_x]=multi_sdr_imperfect(P,h1,h2,rho,beta);
fprintf('multi_im: MA=%d , t=%d, P1=%d, P2=%d Pc=%d Rs=%d \n ',MA_x,t_x, P1_x,P2_x, Pc_x,Rs_x) 

 
fprintf('\n===========================\n\n')

[MA_x,t_x, P1_x,P2_x, Pc_x,Rs_x,fc_x]=multi_kkt_perfect(P,h1,h2,rho);
fc_x;
fprintf('multi: MA=%d , t=%d, P1=%d, P2=%d Pc=%d Rs=%d \n ',MA_x,t_x, P1_x,P2_x, Pc_x,Rs_x)
[MA_x,t_x, P1_x,P2_x, Pc_x,Rs_x,fc_x]=multi_kkt_imperfect(P,h1,h2,rho,beta);
fc_x;
fprintf('multi_im: MA=%d , t=%d, P1=%d, P2=%d Pc=%d Rs=%d \n ',MA_x,t_x, P1_x,P2_x, Pc_x,Rs_x) 


fprintf('\n===========================\n\n')


[MA_x,t_x, P1_x,P2_x, Pc_x,Rs_x,fc_x]=OMA_perfect(P,h1,h2,rho);

fprintf('OMA: MA=%d , t=%d, P1=%d, P2=%d Pc=%d Rs=%d \n ',MA_x,t_x, P1_x,P2_x, Pc_x,Rs_x)
[MA_x,t_x, P1_x,P2_x, Pc_x,Rs_x,fc_x]=OMA_imperfect(P,h1,h2,rho,beta);

fprintf('OMA_im: MA=%d , t=%d, P1=%d, P2=%d Pc=%d Rs=%d \n ',MA_x,t_x, P1_x,P2_x, Pc_x,Rs_x) 

