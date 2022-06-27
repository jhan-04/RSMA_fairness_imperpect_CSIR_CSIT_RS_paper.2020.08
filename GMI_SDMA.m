function  [p_nc,RRR_sdma]=GMI_SDMA(Nt,K,H_h,Pt,N0,sigma_e)


[A,B,C,D]=cal_ABCD(Nt,K,Pt,N0,H_h,sigma_e);
[X_nc,result_nc]=cal_X_no_common(Nt,K,Pt,A,B,C,D);
%result_nc_set(1:length(result_nc))=result_nc;
%[GMI_X_nc]=cal_GMI_withX(K,A,B,C,D,X_nc);
[p_nc]=find_p(K,Pt,A,B,C,D,X_nc);norm(p_nc)^2;
 [GMI,GMI_c,GMI_p]=cal_GMI(K,A,B,C,D,p_nc);


RRR_sdma= [GMI,GMI_c,GMI_p]';
end