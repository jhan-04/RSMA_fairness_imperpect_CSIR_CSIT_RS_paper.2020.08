function  [X,CC,result_set]=cal_X(Nt,K,Pt,A,B,C,D)



%initialization
%Xs(:,:,1)=eye(Nt*(M+1),Nt*(M+1));
%a0(1)=0;
b0(:,1)=(rand(K,1));%randn(M,1)*0+5;%
% %c0(1)=0;
d0(:,1)=(rand(K,1));%randn(M,1)*0+5;%
%bb=(rand(Nt*(M+1),1));%randn(M,1)*0+5;%
%c0(1)=0;
%dd=(rand(Nt*(M+1),1));%randn(M,1)*0+5;%


s=1;
result(1)=0;
GMI(1)=-10;

stop=0;
ee=0.001;

%repeat
while ~stop
    
    s=s+1;
    if b0(1,s-1)~=b0(1,s-1)
        b0(:,s-1)=(rand(K,1));
        d0(:,s-1)=(rand(K,1));
    end
    
     cvx_begin quiet
    %cvx_begin
    variable X(Nt*(K+1),Nt*(K+1)) complex semidefinite
    variable Rt(1)%Lc(1)
    variable a(K,1)%common rate
    variable b(K,1)%common rate
    variable c(K,1)%private rate
    variable d(K,1)%private rate
    variable CC(K,1)%%common part
    %xx=Lc+sum(c)-sum(d);
    maximize(Rt)
    subject to
    for i=1:K
        %a(i)-b(i)>=Lc;
        trace(A(:,:,i)*X)>=exp(a(i));
        trace(C(:,:,i)*X)>=exp(c(i));
        trace(B(:,:,i)*X)<=exp(b0(i,s-1))*(b(i)-b0(i,s-1)+1);
        trace(D(:,:,i)*X)<=exp(d0(i,s-1))*(d(i)-d0(i,s-1)+1);
        
        
        %%add for max min
        a(i)-b(i)>=sum(CC);
        c(i)-d(i)+CC(i)>=Rt;
        
        
    end
    %%add for max min
    CC>=0;
    X== hermitian_semidefinite(Nt*(K+1));
    
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
    [p_max]=find_p(K,Pt,A,B,C,D,X,CC);%norm(p_max(:,i,j))^2;
    [GMI(s),GMI_c,GMI_p]=cal_GMI_min(K,A,B,C,D,p_max,CC);
%     if abs(GMI(s)-GMI(s-1))<ee
%         stop=1;
%     end
    
end

% fprintf('iteration with SDR and cccp  %d in Rs \n',s-2)
% %GMI
% %result/log(2)
%[GMI(s)]=cal_GMI_withX(K,A,B,C,D,X);

result_final=result(s);
result_set=result.';

end