function  [RRR_RS,RRR_sdma,RRR_zf,RRR_mrt,RRR_noma]=CAL_SR_poroposed(Nt,K,H_h,Pt,N0,sigma_e)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Optimization GMI
%%%%%%%%%%%Rate splitting%%%%%%%%%%%%%%
[P_RS,RRR_RS]=GMI_RS(Nt,K,H_h,Pt,N0,sigma_e);%proposed


%%%%%%%%%%%%no_common_part(=SDMA)%%%%%%%%%%%%%%
[P_sdma,RRR_sdma]=GMI_SDMA(Nt,K,H_h,Pt,N0,sigma_e);

%%%%%%%%%Zero_Forcing%%%%%%%%%%
[P_zf,RRR_zf]=GMI_RS_ZF(Nt,K,H_h,Pt,N0,sigma_e);

%%%%%%%%%MRT%%%%%%%%%%
[P_mrt,RRR_mrt]=GMI_RS_MRT(Nt,K,H_h,Pt,N0,sigma_e);

%%%%%%NOMA_two user
[P_noma,RRR_noma]=GMI_NOMA_two_user(Nt,K,H_h,Pt,N0,sigma_e);



end