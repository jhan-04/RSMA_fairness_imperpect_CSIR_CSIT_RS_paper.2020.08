
function  [A,B,C,D]=cal_ABCD(Nt,K,Pt,N0,h,sigma_e)




%%%%A_k
for k=1:K
    H(:,:,k)=h(:,k)*h(:,k)';
    for j=1:K+1
        a(Nt*(j-1)+1:j*Nt,Nt*(j-1)+1:j*Nt,k)=H(:,:,k);
    end
    %A(:,:,k)=[]
    A(:,:,k)=a(:,:,k)+sigma_e(k)^2*eye(Nt*(K+1))+(N0/Pt)*eye(Nt*(K+1));
end

%%%%B_k
for k=1:K
    %H(:,:,k)=h(:,k)*h(:,k)';
    for j=1:K
        b(Nt*(j-1)+1:j*Nt,Nt*(j-1)+1:j*Nt,k)=H(:,:,k);
    end
    b(Nt*(K+1),Nt*(K+1),k)=0;
    %A(:,:,k)=[]
    B(:,:,k)=b(:,:,k)+sigma_e(k)^2*eye(Nt*(K+1))+(N0/Pt)*eye(Nt*(K+1));
end

%%%%C_k
C=B;

%%%%D_k
for k=1:K
    %H(:,:,k)=h(:,k)*h(:,k)';
    for j=1:K
        
        d(Nt*(j-1)+1:j*Nt,Nt*(j-1)+1:j*Nt,k)=H(:,:,k);
        
        if j==k
            d(Nt*(j-1)+1:j*Nt,Nt*(j-1)+1:j*Nt,k)=0;
        end
        
    end
    d(Nt*(K+1),Nt*(K+1),k)=0;
    %A(:,:,k)=[]
    D(:,:,k)=d(:,:,k)+sigma_e(k)^2*eye(Nt*(K+1))+(N0/Pt)*eye(Nt*(K+1));
end

A(:,:,:)=round(A(:,:,:),13);
B(:,:,:)=round(B(:,:,:),13);
C(:,:,:)=round(C(:,:,:),13);
D(:,:,:)=round(D(:,:,:),13);



end
