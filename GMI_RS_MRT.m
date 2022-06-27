function  [p_mrt,RRR_mrt]=GMI_RS_MRT(Nt,K,H_h,Pt,N0,sigma_e)


[A,B,C,D]=cal_ABCD(Nt,K,Pt,N0,H_h,sigma_e);
[Xc,PK,result_mrt]=cal_X_MRT(Nt,K,Pt,H_h,N0,sigma_e);

for k=1:K
    mag(k)=norm(H_h(:,k));
end
pk=reshape(H_h.*(sqrt(PK)./mag')',[],1);
[p_mrt]=find_p_ZF(K,Pt,A,B,C,D,Xc,pk);norm(p_mrt)^2;
[GMI,GMI_c,GMI_p]=cal_GMI(K,A,B,C,D,p_mrt);
RRR_mrt=[GMI,GMI_c,GMI_p]';

end

