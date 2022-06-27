function  [SR,Rs_c,Rs_k]=Cal_Rate(pc,pk,H,K,N0)



for k=1:K
    T_ck(k)=abs(H(:,k)'*pc)^2+sum(abs(H(:,k)'*pk).^2)+N0;
    T_k(k)=sum(abs(H(:,k)'*pk).^2)+N0;
    S_ck(k)=abs(H(:,k)'*pc)^2;
    S_k(k)=abs(H(:,k)'*pk(:,k)).^2;
    I_ck(k)=T_k(k);
    I_k(k)=I_ck(k)-S_k(k);
    Rs_ck(k)=log2(1+S_ck(k)/I_ck(k));
    Rs_k(k)=log2(1+S_k(k)/I_k(k));
    

    
end
%%%%
Rs_c=min(Rs_ck);
SR=Rs_c+sum(Rs_k);


end