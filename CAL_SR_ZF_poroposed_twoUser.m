function  [RRR_zf_i,RRR_zf_p,RRR_zf_ni,RRR_zf_np,RRR_zf_si,RRR_zf_sp,RRR_zf_mi,RRR_zf_mp,RRR_mrt1,RRR_mrt2]=CAL_SR_ZF_poroposed_twoUser(H_h,Pt,beta)

i=1;
j=1;
if norm(H_h(:,1))>=norm(H_h(:,2))
    h1(:,i)=H_h(:,1); h2(:,i)=H_h(:,2);
else
    h1(:,i)=H_h(:,2); h2(:,i)=H_h(:,1);
end

rho=1-abs(h1(:,i)'/norm(h1(:,i))*h2(:,i)/norm(h2(:,i)))^2;

%%%%%%%%%RS_NOMA
[MA_Ri(i,j),t_Ri(i,j), P1_Ri(i,j),P2_Ri(i,j), Pc_Ri(i,j),Rs_Ri(i,j),fc_Ri(:,i,j),t_Ni(i,j),Rs_Ni(i,j),fc_Ni(:,i,j)]=RSNOMA_sdr_imperfect(Pt,h1(:,i),h2(:,i),rho,beta);
[MA_Rp(i,j),t_Rp(i,j), P1_Rp(i,j),P2_Rp(i,j), Pc_Rp(i,j),Rs_Rp(i,j),fc_Rp(:,i,j),t_Np(i,j),Rs_Np(i,j),fc_Np(:,i,j)]=RSNOMA_sdr_imperfect(Pt,h1(:,i),h2(:,i),rho,0);
%%%%%%%%%%%%%%%%%%%%%%%%%SDMA
[MA_Si(i,j),t_Si(i,j), P1_Si(i,j),P2_Si(i,j), Pc_Si(i,j),Rs_Si(i,j),fc_Si(:,i,j)]=SDMA_imperfect(Pt,h1(:,i),h2(:,i),rho,beta);
[MA_Sp(i,j),t_Sp(i,j), P1_Sp(i,j),P2_Sp(i,j), Pc_Sp(i,j),Rs_Sp(i,j),fc_Sp(:,i,j)]=SDMA_perfect(Pt,h1(:,i),h2(:,i),rho);
%%%%%%%MRT
[MA_Oim1(i,j),t_Oim1(i,j), P1_Oim1(i,j),P2_Oim1(i,j), Pc_Oim1(i,j),Rs_Oim1(i,j),fc_Oim1(:,i,j)]=OMA_MRT_one_imperfect(Pt,h1(:,i),h2(:,i),rho,beta);
[Rs_x,Rs_x1,Rs_x2]=OMA_MRT_imperfect(Pt,h1(:,i),h2(:,i),rho,beta);%two user
%%%%%%multicasting
[MA_mi(i,j),t_mi(i,j), P1_mi(i,j),P2_mi(i,j), Pc_mi(i,j),Rs_mi(i,j),fc_mi(:,i,j)]=multi_sdr_imperfect(Pt,h1(:,i),h2(:,i),rho,beta);
[MA_mp(i,j),t_mp(i,j), P1_mp(i,j),P2_mp(i,j), Pc_mp(i,j),Rs_mp(i,j),fc_mp(:,i,j)]=multi_sdr_perfect(Pt,h1(:,i),h2(:,i),rho);

%%%%%%%%%%%%%%%%%%=====caculate GMI=======
%%%%%RS NOMA
[Rs_GMI,Rc_x,R1_x,R2_x]=cal_GMI_1(Pt,t_Ri(i,j), P1_Ri(i,j),P2_Ri(i,j), Pc_Ri(i,j),fc_Ri(:,i,j),h1(:,i),h2(:,i),rho,beta);
RRR_zf_i=[Rs_GMI,Rc_x,R1_x,R2_x]';
[Rs_GMI,Rc_x,R1_x,R2_x]=cal_GMI_1(Pt,t_Rp(i,j), P1_Rp(i,j),P2_Rp(i,j), Pc_Rp(i,j),fc_Rp(:,i,j),h1(:,i),h2(:,i),rho,beta);
RRR_zf_p=[Rs_GMI,Rc_x,R1_x,R2_x]';

[Rs_GMI,Rc_x,R1_x,R2_x]=cal_GMI_1(Pt,t_Ni(i,j), t_Ni(i,j)*Pt,0, (1-t_Ni(i,j))*Pt,fc_Ni(:,i,j),h1(:,i),h2(:,i),rho,beta);
RRR_zf_ni=[Rs_GMI,Rc_x,R1_x,R2_x]';
[Rs_GMI,Rc_x,R1_x,R2_x]=cal_GMI_1(Pt,t_Np(i,j), t_Np(i,j)*Pt,0, (1-t_Np(i,j))*Pt,fc_Np(:,i,j),h1(:,i),h2(:,i),rho,beta);
RRR_zf_np=[Rs_GMI,Rc_x,R1_x,R2_x]';
%%%%%%%%%%%%%%%%%%%%%%%%%SDMA
[Rs_GMI,Rc_x,R1_x,R2_x]=cal_GMI_1(Pt,t_Si(i,j), P1_Si(i,j),P2_Si(i,j), Pc_Si(i,j),fc_Si(:,i,j),h1(:,i),h2(:,i),rho,beta);
RRR_zf_si=[Rs_GMI,Rc_x,R1_x,R2_x]';
[Rs_GMI,Rc_x,R1_x,R2_x]=cal_GMI_1(Pt,t_Sp(i,j), P1_Sp(i,j),P2_Sp(i,j), Pc_Sp(i,j),fc_Sp(:,i,j),h1(:,i),h2(:,i),rho,beta);
RRR_zf_sp=[Rs_GMI,Rc_x,R1_x,R2_x]';
%%%%%%%%%%%%%%%multi
[Rs_GMI,Rc_x,R1_x,R2_x]=cal_GMI_1(Pt,t_mi(i,j), P1_mi(i,j),P2_mi(i,j), Pc_mi(i,j),fc_mi(:,i,j),h1(:,i),h2(:,i),rho,beta);
RRR_zf_mi=[Rs_GMI,Rc_x,R1_x,R2_x]';
[Rs_GMI,Rc_x,R1_x,R2_x]=cal_GMI_1(Pt,t_mp(i,j), P1_mp(i,j),P2_mp(i,j), Pc_mp(i,j),fc_mp(:,i,j),h1(:,i),h2(:,i),rho,beta);
RRR_zf_mp=[Rs_GMI,Rc_x,R1_x,R2_x]';
%%%%MRT
RRR_mrt1=[Rs_Oim1(i,j),0,Rs_Oim1(i,j),0]';
RRR_mrt2=[Rs_x,0,Rs_x1,Rs_x2]';


end