function [GMI_min,GMI_c,GMI_p]=cal_GMI_min(K,A,B,C,D,p,CC)


for k=1:K
    GMI_p(k)=real(log2((p'*C(:,:,k)*p)/(p'*D(:,:,k)*p)));
    GMI_cc(k)=real(log2((p'*A(:,:,k)*p)/(p'*B(:,:,k)*p)));
end
GMI_c=min(GMI_cc);
GMI=sum(GMI_p)+min(GMI_cc);

GMI_min=min(GMI_p+GMI_c.*CC'./sum(CC));

end