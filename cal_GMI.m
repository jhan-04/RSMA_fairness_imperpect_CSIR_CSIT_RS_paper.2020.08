function [GMI,GMI_c,GMI_p]=cal_GMI(K,A,B,C,D,p)


for k=1:K
    GMI_p(k)=log2((p'*C(:,:,k)*p)/(p'*D(:,:,k)*p));
    GMI_cc(k)=log2((p'*A(:,:,k)*p)/(p'*B(:,:,k)*p));
end
GMI_c=min(GMI_cc);
GMI=sum(GMI_p)+min(GMI_cc);

end