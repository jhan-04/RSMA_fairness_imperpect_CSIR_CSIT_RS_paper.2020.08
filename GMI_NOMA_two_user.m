function  [p_N,RRR_noma]=GMI_NOMA_two_user(Nt,K,H_h,Pt,N0,sigma_e)


        
        [A,B,C,D]=cal_ABCD(Nt,K,Pt,N0,H_h,sigma_e);
        [X_N,result_set_N]=cal_X_NOMA3(Nt,K,H_h,A,B,C,D);
        [p_N]=find_p_NOMA(K,Pt,A,B,C,D,X_N,H_h);
        [GMI,GMI_c,GMI_p]=cal_GMI(K,A,B,C,D,p_N);

RRR_noma=[GMI,GMI_c,GMI_p];

end
