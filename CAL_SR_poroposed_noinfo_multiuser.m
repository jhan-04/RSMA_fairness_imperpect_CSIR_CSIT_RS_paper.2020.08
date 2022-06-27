function  [RRR_G,RRR_sdma,RRR_zf,RRR_mrt]=CAL_SR_poroposed_noinfo_multiuser(Nt,K,H_h,H_m,Pt,M,N0,N02,sigma_e)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Optimization GMI
%%%%%%%%%%%Rate splitting%%%%%%%%%%%%%%
[P_G,pc_G,pk_G,GMI,SR]=GMI_RS(Nt,K,H_h,H_m,Pt,M,N02,sigma_e);%proposed
[SR,Rs_c_G,Rs_k_G]=Cal_Rate(pc_G,pk_G,H_m,K,N0);
RRR_G=[SR,Rs_c_G,Rs_k_G]';

%%%%%%%%%%%%no_common_part(=SDMA)%%%%%%%%%%%%%%
[P_sdma,pc_sdma,pk_sdma,GMI_sdma,SR_sdma]=GMI_SDMA(Nt,K,H_h,H_m,Pt,M,N02,sigma_e);
[SR_sdma,Rs_c_sdma,Rs_k_sdma]=Cal_Rate(pc_sdma,pk_sdma,H_m,K,N0);
RRR_sdma=[SR_sdma,Rs_c_sdma,Rs_k_sdma]';
%%%%%%%%%Zero_Forcing%%%%%%%%%%
[P_zf,pc_zf,pk_zf,GMI_zf,SR_zf]=GMI_RS_ZF(Nt,K,H_h,H_m,Pt,M,N02,sigma_e);
[SR_zf,Rs_c_zf,Rs_k_zf]=Cal_Rate(pc_zf,pk_zf,H_m,K,N0);
RRR_zf=[SR_zf,Rs_c_zf,Rs_k_zf]';
%%%%%%%%%MRT%%%%%%%%%%
[P_mrt,pc_mrt,pk_mrt,GMI_mrt,SR_mrt]=GMI_RS_MRT(Nt,K,H_h,H_m,Pt,M,N02,sigma_e);
[SR_mrt,Rs_c_mrt,Rs_k_mrt]=Cal_Rate(pc_mrt,pk_mrt,H_m,K,N0);
RRR_mrt=[SR_mrt,Rs_c_mrt,Rs_k_mrt]';



end