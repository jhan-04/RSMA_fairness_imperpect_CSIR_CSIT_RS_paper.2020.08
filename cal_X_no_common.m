function  [X_nc,result_set]=cal_X_no_common(Nt,K,Pt,A,B,C,D)




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%initialization
%Xs(:,:,1)=eye(Nt*(M+1),Nt*(M+1));
%a0(1)=0;
%b0(:,1)=abs(randn(M,1))*3;
%c0(1)=0;
d0(:,1)=rand(K,1);

s=1;
result(1)=0;

stop=0;
ee=0.001;
for k=1:K
    CC(:,:,k)=C(1:Nt*K,1:Nt*K,k);
    DD(:,:,k)=D(1:Nt*K,1:Nt*K,k);
end

%repeat
while ~stop
    
    s=s+1;
    if d0(1,s-1)~=d0(1,s-1)
        d0(:,s-1)=(rand(K,1));
    end
    
    cvx_begin quiet
    %cvx_begin
    %variable X(Nt*(K+1),Nt*(K+1)) hermitian
    variable X_nc(Nt*(K),Nt*(K)) hermitian
    variable Rt(1)%Lc(1) nonnegative
    %     variable a(M,1) nonnegative
    %     variable b(M,1) nonnegative
    variable c(K,1) nonnegative
    variable d(K,1) nonnegative
    %xx=sum(c)-sum(d);
    maximize(Rt)
    
    subject to
    for i=1:K
        %a(i)-b(i)>=Lc
        %trace(A(:,:,i)*X)>=exp(a(i))
        trace(CC(:,:,i)*X_nc)>=exp(c(i))
        %linearize
        %trace(B(:,:,i)*X)<=exp(b0(i,s-1))*(b(i)-b0(i,s-1)+1)
        trace(DD(:,:,i)*X_nc)<=exp(d0(i,s-1))*(d(i)-d0(i,s-1)+1)
        
         %%add for max min
        c(i)-d(i)>=Rt;
    end
    %trace(X_nc)<=Pt;
    X_nc==hermitian_semidefinite(Nt*(K));
    
    cvx_end
    
    %     b0(:,s)=b;
    d0(:,s)=d;
    
    
    % until
    if norm(d0(:,s)-d0(:,s-1))^2<ee %norm(result(s)-result(s-1))^2<0.0001
        stop=1;
    end
    
    
    result(s)=Rt;
    if result(s-1)>result(s)&&s>2
        fprintf('NOT monotonically increasing %d in No- Rs \n',result(s-1)-result(s))
        result
    end
    X_nc1=X_nc;
    X_nc(Nt*(K+1),Nt*(K+1))=0;
    
    
    
    for k=1:K
        GMI_p(k)=log2(trace(C(:,:,k)*X_nc))-log2((trace(D(:,:,k)*X_nc)));
        %GMI_c(k)=log2(trace(A(:,:,k)*X_nc)/(trace(B(:,:,k)*X_nc)));
    end
    
    GMI(s)=sum(GMI_p);
    
    
    
end

GMI;
result/log(2);
result_final=result(s);
result_set=result.';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



end