function  [RRR_RSe,RRR_sdmae,RRR_zfe,RRR_mrte,RRR_nomae]=CAL_SR_poroposed_noinfo(Nt,K,H_h,Pt,N0,sigma_e)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Optimization GMI

        [A,B,C,D]=cal_ABCD(Nt,K,Pt,N0,H_h,sigma_e);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Optimization GMI
%%%%%%%%%%%Rate splitting%%%%%%%%%%%%%%
[P_RS,RRR_RS]=GMI_RS(Nt,K,H_h,Pt,N0,sigma_e*0);%proposed
[GMI,GMI_c,GMI_p]=cal_GMI(K,A,B,C,D,P_RS);
RRR_RSe=[GMI,GMI_c,GMI_p]';
%%%%%%%%%%%%no_common_part(=SDMA)%%%%%%%%%%%%%%
[P_sdma,RRR_sdma]=GMI_SDMA(Nt,K,H_h,Pt,N0,sigma_e*0);
[GMI,GMI_c,GMI_p]=cal_GMI(K,A,B,C,D,P_sdma);
RRR_sdmae=[GMI,GMI_c,GMI_p]';
%%%%%%%%%Zero_Forcing%%%%%%%%%%
[P_zf,RRR_zf]=GMI_RS_ZF(Nt,K,H_h,Pt,N0,sigma_e*0);
[GMI,GMI_c,GMI_p]=cal_GMI(K,A,B,C,D,P_zf);
RRR_zfe=[GMI,GMI_c,GMI_p]';
%%%%%%%%%MRT%%%%%%%%%%
[P_mrt,RRR_mrt]=GMI_RS_MRT(Nt,K,H_h,Pt,N0,sigma_e*0);
[GMI,GMI_c,GMI_p]=cal_GMI(K,A,B,C,D,P_mrt);
RRR_mrte=[GMI,GMI_c,GMI_p]';
%%%%%%NOMA_two user
[P_noma,RRR_noma]=GMI_NOMA_two_user(Nt,K,H_h,Pt,N0,sigma_e*0);
[GMI,GMI_c,GMI_p]=cal_GMI(K,A,B,C,D,P_noma);
RRR_nomae=[GMI,GMI_c,GMI_p]';
end