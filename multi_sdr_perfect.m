function [MA_x,t_x, P1_x,P2_x, Pc_x,Rs_x,fc_x]=multi_sdr_perfect(P,h1,h2,rho)


MA_x=4;
P2_x=0;
P1_x=0;
Pc_x=P;
t_x=0;
Nt=length(h1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A=h1*h1';
B=h2*h2';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cvx_begin quiet
variable F(Nt,Nt) hermitian
variable x(1)
minimize(-x)
subject to
trace(A*F)>=x;
trace(B*F)>=x;
trace(F)==1;
F==hermitian_semidefinite(Nt);
cvx_end
%%%%%%%%%%%%%%%%%%%%%%%%%%

trace(B*F);

[U,W,Z] = svds(F);%V=U*W*Z'
obj_max=-100;
n=length(W);
for m=1:10000
    r=sqrt(1/2)*(randn(n,1)+1i*randn(n,1));%sqrt(var/2)*(randn(1,N)+1i*randn(1,N))
    v_=U*sqrt(W)*r;
    v_=v_/sqrt(v_'*v_);
    obj=min(v_'*A*v_,v_'*B*v_);
    if obj>=obj_max
        obj_max=obj;
        v_max=v_;
    end
end
fc_x=v_max;


Rs_x=min(log2(1+abs(h1'*fc_x)^2*P),log2(1+abs(h2'*fc_x)^2*P));


%Rc1_mul_m=log2(1+P*abs(h1'*f_c)^2);
%Rc2_mul_m=log2(1+P*abs(h2'*f_c)^2);


end