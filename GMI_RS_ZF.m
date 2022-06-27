function  [p_zf,RRR_ZF]=GMI_RS_ZF(Nt,K,H_h,Pt,N0,sigma_e)


[A,B,C,D]=cal_ABCD(Nt,K,Pt,N0,H_h,sigma_e);
[Xc,PK,result]=cal_X_ZF(Nt,K,Pt,H_h,N0,sigma_e);
%result_set(1:length(result))=result;
pk_zf=H_h*inv(H_h'*H_h);
for k=1:K
    mag(k)=norm(pk_zf(:,k));
end
pk=reshape(pk_zf.*(sqrt(PK)./mag')',[],1);
[p_zf]=find_p_ZF(K,Pt,A,B,C,D,Xc,pk);norm(p_zf)^2;
[GMI,GMI_c,GMI_p]=cal_GMI(K,A,B,C,D,p_zf);
RRR_ZF=[GMI,GMI_c,GMI_p]';
end
