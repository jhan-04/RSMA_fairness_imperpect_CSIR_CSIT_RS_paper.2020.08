function  [Xc_final,P_final,result_final]=cal_X_MRT(Nt,M,Pt,h,N0,sigma_e)



%h=Nt X M matrix
% NtXM matrix

for k=1:M
    for i=1:M
        A(i,k)=abs(h(:,k)'*h(:,i)/norm(h(:,i)))^2;
    end
end


test=1;
for t=1:test
    
    
    %initialization
    %Xs(:,:,1)=eye(Nt*(M+1),Nt*(M+1));
    %a0(1)=0;
    b0(:,1)=log((rand(M,1))*(Pt/M)+(sigma_e.^2*Pt+N0));%randn(M,1)*0+5;%
    %b0(:,1)=(rand(M,1))*log(Pt/M)+log(sigma_e.^2*Pt+N0);%randn(M,1)*0+5;%
    %c0(1)=0;
    %d0(:,1)=(rand(M,1))*log(Pt/M)+log(sigma_e.^2*Pt+N0);%randn(M,1)*0+5;%
    d0(:,1)=log((rand(M,1))*(Pt/M)+(sigma_e.^2*Pt+N0));%randn(M,1)*0+5;%
%     b0(:,1)
%     d0(:,1)
    
    s=1;
    result(1)=0;
    
    stop=0;
    ee=0.001;
    
    
    %repeat
    while ~stop
        
        s=s+1;
        %b0(:,s-1)
        if b0(1,s-1)~=b0(1,s-1)
            b0(:,s-1)=log((rand(M,1))*(Pt/M)+(sigma_e.^2*Pt+N0));
            d0(:,s-1)=log((rand(M,1))*(Pt/M)+(sigma_e.^2*Pt+N0));
        end
        
        cvx_begin quiet
        %cvx_begin
        variable Xc(Nt,Nt) hermitian
        variable Lc(1) nonnegative
        variable a(M,1) nonnegative
        variable b(M,1) nonnegative
        variable c(M,1) nonnegative
        variable d(M,1) nonnegative
        variable P(M,1) nonnegative
        xx=Lc+sum(c)-sum(d);
        maximize(xx)
        subject to
        
        for i=1:M
            
            a(i)-b(i)>=Lc;
            H=(h(:,i)*h(:,i)'+(h(:,i)*h(:,i)')')/2;
            
            %AA(i)=aa(i)*P(i);
            trace(H*Xc)+ A(:,i)'*P+sigma_e(i)^2*Pt+N0>=exp(a(i));
            A(:,i)'*P+sigma_e(i)^2*sum(P)+N0>=exp(c(i));
            
            A(:,i)'*P+sigma_e(i)^2*Pt+N0<=exp(b0(i,s-1))*(b(i)-b0(i,s-1)+1);
            A(:,i)'*P-norm(h(:,i))^2*P(i)+sigma_e(i)^2*sum(P)+N0<=exp(d0(i,s-1))*(d(i)-d0(i,s-1)+1);
            
            
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
    
    
    test_result(t)=result(s);
    test_Xc(:,:,t)=Xc;
    test_P(:,t)=P;
end
% %GMI
% %result/log(2)
result_final=max(test_result);
tt=find(max(test_result)==test_result);
Xc_final=test_Xc(:,:,tt);
P_final=test_P(:,tt);
%result_final=result(s);
result_set=result.';
cd=sum(c)-sum(d);

end