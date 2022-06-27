function  [Xc,P,result_set]=cal_X_ZF(Nt,K,Pt,h,N0,sigma_e)



%h=Nt X M matrix
p=h*inv(h'*h);% NtXM matrix

for k=1:K
    A(:,k)=zeros(K,1);
    A(k,k)=abs(h(:,k)'*p(:,k)/norm(p(:,k)))^2;
    aa(k)=norm(h(:,k)'*p(:,k)/norm(p(:,k)))^2;
end

%initialization
%Xs(:,:,1)=eye(Nt*(M+1),Nt*(M+1));
%a0(1)=0;
b0(:,1)=log((rand(K,1)).*(Pt/K)+(sigma_e.^2*Pt+N0));%randn(K,1)*0+5;%
%c0(1)=0;
d0(:,1)=log((rand(K,1)).*(Pt/K)+(sigma_e.^2*Pt+N0));%randn(K,1)*0+5;%

s=1;
result(1)=0;

stop=0;
ee=0.001;

%repeat
while ~stop
    
    s=s+1;
     %b0(:,s-1)
    if b0(1,s-1)~=b0(1,s-1)
        b0(:,s-1)=log((rand(K,1))*(Pt/K)+(sigma_e.^2*Pt+N0));
        d0(:,s-1)=log((rand(K,1))*(Pt/K)+(sigma_e.^2*Pt+N0));
    end
   
    cvx_begin quiet
    %cvx_begin
    variable Xc(Nt,Nt) hermitian
    variable Lc(1) nonnegative
    variable a(K,1) nonnegative
    variable b(K,1) nonnegative
    variable c(K,1) nonnegative
    variable d(K,1) nonnegative
    variable P(K,1) nonnegative
    xx=Lc+sum(c)-sum(d);
    maximize(xx)
    subject to
    
    for i=1:K
        
        a(i)-b(i)>=Lc;
        H=round(h(:,i)*h(:,i)',10);
        
        %AA(i)=aa(i)*P(i);
        trace(H*Xc)+aa(i)*P(i)+sigma_e(i)^2*Pt+N0>=exp(a(i));
        aa(i)*P(i)+sigma_e(i)^2*sum(P)+N0>=exp(c(i));
        
        aa(i)*P(i)+sigma_e(i)^2*Pt+N0<=exp(b0(i,s-1))*(b(i)-b0(i,s-1)+1);
        sigma_e(i)^2*sum(P)+N0<=exp(d0(i,s-1))*(d(i)-d0(i,s-1)+1);
        
        
    end
    
    
    trace(Xc)+sum(P)==Pt;
    Xc==hermitian_semidefinite(Nt);
    
    cvx_end
    
    
    %zzz=[a b c d]
    
    b0(:,s)=b;
    d0(:,s)=d;
    result(s)=Lc+sum(c)-sum(d);
    result_Lc(s)=Lc;
    result_c(s)=sum(c);
    result_d(s)=sum(d);
    
    % until
    if norm(b0(:,s)-b0(:,s-1))^2+norm(d0(:,s)-d0(:,s-1))^2<ee%norm(result(s)-result(s-1))^2<0.0001
        stop=1;
    end
    
    if result(s-1)>result(s)&&s>2
        fprintf('NOT monotonically increasing %d in Rs_ZF \n',result(s-1)-result(s))
        result
    end
    %
    
end

% %GMI
% %result/log(2)
result_final=result(s);
result_set=result.';
cd=sum(c)-sum(d);

end