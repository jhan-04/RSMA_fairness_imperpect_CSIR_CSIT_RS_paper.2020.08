function  [p_opt,result_cccp,GMI_cccp]=GMI_only_cccp_ver2(Nt,K,H_h,Pt,N0,sigma_e,num_interation)

[A,B,C,D]=cal_ABCD(Nt,K,Pt,N0,H_h,sigma_e);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%only CCCP
%initialization

b0(:,1)=(rand(K,1));%randn(M,1)*0+5;%
d0(:,1)=(rand(K,1));%randn(M,1)*0+5;%

p0(:,1)=(rand(Nt*(K+1),1)).*exp(1i.*2.*pi.*rand(Nt*(K+1),1));%randn(M,1)*0+5;%
p0(:,1)=p0(:,1)/norm(p0(:,1));

GMI_cccp(1)=0;
s=1;
result_cccp(1)=0;

stop=0;
ee=0.001;

%repeat
while ~stop
    
    s=s+1;
    
    cvx_begin quiet
    %cvx_begin
    variable p(Nt*(K+1),1) complex
    variable p_comp(Nt*(K+1),1)
    variable p_real(Nt*(K+1),1)
    variable Lc(1)
    variable a(K,1)
    variable b(K,1)
    variable c(K,1)
    variable d(K,1)
    xx=Lc+sum(c)-sum(d);
    maximize(xx)
    subject to
    for i=1:K
        a(i)-b(i)>=Lc;
        %   p==1i*p_comp+p_real;
        quad_form(p,B(:,:,i))<=exp(b0(i,s-1))*(b(i)-b0(i,s-1)+1);
        quad_form(p,D(:,:,i))<=exp(d0(i,s-1))*(d(i)-d0(i,s-1)+1);
        [real(2.*p0(:,s-1)'*A(:,:,i)), -imag(2.*p0(:,s-1)'*A(:,:,i))]*[real(p);imag(p)]-real(p0(:,s-1)'*A(:,:,i)*p0(:,s-1))>=exp(a(i));
        [real(2.*p0(:,s-1)'*C(:,:,i)), -imag(2.*p0(:,s-1)'*C(:,:,i))]*[real(p);imag(p)]-real(p0(:,s-1)'*C(:,:,i)*p0(:,s-1))>=exp(c(i));
        %
        [real(2.*p0(:,s-1)'*A(:,:,i)), imag(2.*p0(:,s-1)'*A(:,:,i))]*[imag(p);real(p)]>=0;
        [real(2.*p0(:,s-1)'*C(:,:,i)), imag(2.*p0(:,s-1)'*C(:,:,i))]*[imag(p);real(p)]>=0;
%         [real(p0(:,s-1)'*A(:,:,i)), -imag(p0(:,s-1)'*A(:,:,i))]*[real(p);imag(p)]>=exp(a(i));
%         [real(p0(:,s-1)'*C(:,:,i)), -imag(p0(:,s-1)'*C(:,:,i))]*[real(p);imag(p)]>=exp(c(i));
%         %
%         [real(p0(:,s-1)'*A(:,:,i)), imag(p0(:,s-1)'*A(:,:,i))]*[imag(p);real(p)]>=0;
%         [real(p0(:,s-1)'*C(:,:,i)), imag(p0(:,s-1)'*C(:,:,i))]*[imag(p);real(p)]>=0;
        
    end
         %norm(p)<=Pt;
    cvx_end
    
    b0(:,s)=b; d0(:,s)=d; p0(:,s)=p;
    result_cccp(s)=Lc+sum(c)-sum(d);
    
    % until
    [GMI_cccp(s),GMI_c_cccp,GMI_p_cccp]=cal_GMI(K,A,B,C,D,p);
    if s>num_interation%abs(GMI_cccp(s)-GMI_cccp(s-1))<ee
        stop=1;
    end
    
    if result_cccp(s-1)>result_cccp(s)&&s>2
        fprintf('NOT monotonically increasing %d in Rs \n',result_cccp(s-1)-result_cccp(s))
        result_cccp
    end
    %
    % [GMI(s)]=cal_GMI_withX(K,A,B,C,D,X);
    %p_opt(:,s)=p;
end
fprintf('iteration with only cccp  %d in Rs \n',s-2)
result_cccp(s);
p_opt=p;

end