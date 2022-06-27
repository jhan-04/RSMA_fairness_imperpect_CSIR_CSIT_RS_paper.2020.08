function [GMI]=cal_GMI_withX(M,A,B,C,D,X)


for k=1:M
    GMI_p(k)=log2(trace(C(:,:,k)*X)/(trace(D(:,:,k)*X)));
    GMI_c(k)=log2(trace(A(:,:,k)*X)/(trace(B(:,:,k)*X)));
end

GMI=sum(GMI_p)+min(GMI_c);

end