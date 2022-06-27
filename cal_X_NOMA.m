function [XX,result_set]=cal_X_NOMA(Nt,K,h,A,B,C,D)
%for only two user

if length(h(1,:))~=2
    fprintf('it doesnt work NOMA\n')
end


%initialization
%Xs(:,:,1)=eye(Nt*(M+1),Nt*(M+1));
%a0(1)=0;
b0(:,1)=(rand(K,1));%randn(M,1)*0+5;%
% %c0(1)=0;
d0(:,1)=(rand(1,1));%randn(M,1)*0+5;%
%bb=(rand(Nt*(M+1),1));%randn(M,1)*0+5;%
%c0(1)=0;
%dd=(rand(Nt*(M+1),1));%randn(M,1)*0+5;%
if norm(h(:,1))>norm(h(:,2))
    user1=1;
    user2=3;
else
    user1=3;
    user2=1;
end
for i=1:K
    A([user2,user2+1],:,i)=zeros(2,Nt*(K+1));
    A(:,[user2,user2+1],i)=zeros(Nt*(K+1),2);
    B([user2,user2+1],:,i)=zeros(2,Nt*(K+1));
    B(:,[user2,user2+1],i)=zeros(Nt*(K+1),2);
    C([user2,user2+1],:,i)=zeros(2,Nt*(K+1));
    C(:,[user2,user2+1],i)=zeros(Nt*(K+1),2);
    D([user2,user2+1],:,i)=zeros(2,Nt*(K+1));
    D(:,[user2,user2+1],i)=zeros(Nt*(K+1),2);
    
end
s=1;
result(1)=0;

stop=0;
ee=0.001;

%repeat
while ~stop
    
    s=s+1;
    if b0(1,s-1)~=b0(1,s-1)
        b0(:,s-1)=(rand(K,1));
        d0(:,s-1)=(rand(1,1));
    end
    
    cvx_begin quiet
    %cvx_begin
    variable X(Nt*(K+1),Nt*(K+1)) complex semidefinite
    variable Lc(1)
    variable Rt(1)
    variable a(K,1)
    variable b(K,1)
    variable c(1,1)
    variable d(1,1)
    xx=Lc+(c)-(d);
    maximize(Rt)
    subject to
    trace(C(:,:,1)*X)>=exp(c);
    trace(D(:,:,1)*X)<=exp(d0(1,s-1))*(d-d0(1,s-1)+1);
    for i=1:K
        a(i)-b(i)>=Lc;
        trace(A(:,:,i)*X)>=exp(a(i));
        trace(B(:,:,i)*X)<=exp(b0(i,s-1))*(b(i)-b0(i,s-1)+1);
        trace((round(h(:,K)*h(:,K)',10))*X(1:Nt-1,1:Nt-1))>=trace((round(h(:,K)*h(:,K)',10))'*X(Nt*K+1:Nt*(K+1),Nt*K+1:Nt*(K+1)));
       
    end
       Lc>=Rt;
       c-d>=Rt;
    
    %trace(X)<=Pt;
    X== hermitian_semidefinite(Nt*(K+1));
    %X([3,4],:)==zeros(2,Nt*(K+1));
    X(:,[3,4])==zeros(Nt*(K+1),2);%since X is hermitian, X(:,[3,4])==zeros contraint also satisfy %X([3,4],:)==zeros
    
    cvx_end
    
    
    %zzz=[a b c d]
    
    b0(:,s)=b;
    d0(:,s)=d;
    result(s)=Rt;
    
    % until
    if norm(b0(:,s)-b0(:,s-1))^2+norm(d0(:,s)-d0(:,s-1))^2<ee%norm(result(s)-result(s-1))^2<0.0001
        stop=1;
    end
    
    if result(s-1)>result(s)&&s>2
        fprintf('NOT monotonically increasing %d in Rs \n',result(s-1)-result(s))
        result
    end
    %
    [GMI(s)]=cal_GMI_withX(K,A,B,C,D,X);
    
end
X([3,4],:)=zeros(2,Nt*(K+1));
X(:,[3,4])=zeros(Nt*(K+1),2);
XX=full(X);
% %GMI
% %result/log(2)
XXresult_final=result(s);
result_set=result.';


end