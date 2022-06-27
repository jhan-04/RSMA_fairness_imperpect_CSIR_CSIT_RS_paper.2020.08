function  [p_max,RRR_RS]=GMI_RS(Nt,K,H_h,Pt,N0,sigma_e)


[A,B,C,D]=cal_ABCD(Nt,K,Pt,N0,H_h,sigma_e);
[X,result]=cal_X(Nt,K,Pt,A,B,C,D);
%result_set(1:length(result))=result;
[p_max]=find_p(K,Pt,A,B,C,D,X);%norm(p_max(:,i,j))^2;

[GMI,GMI_c,GMI_p]=cal_GMI(K,A,B,C,D,p_max);

RRR_RS=[GMI,GMI_c,GMI_p]';
end