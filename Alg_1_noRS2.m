function  [pc,pk]=Alg_1_noRS2(Nt,K,H_h,H_m,Pt,M,N0)
%%%


for k=1:K
    pk(:,k)=sqrt(1/K)*(H_h(:,k)/norm(H_h(:,k)));
end
%pk=sqrt(Pt^(-0.6))*(H_h/norm(H_h));
[X,Y,Z]=svd(H_h);
pc=zeros(Nt,1);
pk=sqrt(Pt)*pk;
n=2;
A=[];
A(2)=-10;
A(1)=10;
e=0.001;
while abs((A(n)-A(n-1))/A(n))>=e
    
    for m=1:M
        
        
        for k=1:K
            T_ck(k,m)=abs(H_m(:,k,m)'*pc)^2+sum(abs(H_m(:,k,m)'*pk).^2)+N0;
            T_k(k,m)=sum(abs(H_m(:,k,m)'*pk).^2)+N0;
            S_ck(k,m)=abs(H_m(:,k,m)'*pc)^2;
            S_k(k,m)=abs(H_m(:,k,m)'*pk(:,k))^2;
            I_ck(k,m)=T_k(k,m);
            I_k(k,m)=I_ck(k,m)-S_k(k,m);
            %%%
            g_ck_m(k,m)=pc'*H_m(:,k,m)/T_ck(k,m);
            g_k_m(k,m)=pk(:,k)'*H_m(:,k,m)/T_k(k,m);
            u_ck_m(k,m)=T_ck(k,m)/I_ck(k,m);
            u_k_m(k,m)=T_k(k,m)/I_k(k,m);
            %%%
            t_ck_m(k,m)=u_ck_m(k,m)*abs(g_ck_m(k,m))^2;
            t_k_m(k,m)=u_k_m(k,m)*abs(g_k_m(k,m))^2;
            psi_ck_m(:,:,k,m)=t_ck_m(k,m)*H_m(:,k,m)*H_m(:,k,m)';
            psi_k_m(:,:,k,m)=t_k_m(k,m)*H_m(:,k,m)*H_m(:,k,m)';
            f_ck_m(:,k,m)=u_ck_m(k,m)*H_m(:,k,m)*g_ck_m(k,m)';
            f_k_m(:,k,m)=u_k_m(k,m)*H_m(:,k,m)*g_k_m(k,m)';
            v_ck_m(k,m)=log2(u_ck_m(k,m));
            v_k_m(k,m)=log2(u_k_m(k,m));
        end
    end
    %averaging=bar{variable in paper}
    for k=1:K
        t_ck(k)=mean(t_ck_m(k,:));
        t_k(k)=mean(t_k_m(k,:));
        psi_ck(:,:,k)=(mean(psi_ck_m(:,:,k,:),4)+mean(psi_ck_m(:,:,k,:),4)')/2;
        psi_k(:,:,k)=(mean(psi_k_m(:,:,k,:),4)+mean(psi_k_m(:,:,k,:),4)')/2;
        f_ck(:,k)=mean(f_ck_m(:,k,:),3);
        f_k(:,k)=mean(f_k_m(:,k,:),3);
        v_ck(k)=mean(v_ck_m(k,:));
        v_k(k)=mean(v_k_m(k,:));
        %%%
        u_ck(k)=mean(u_ck_m(k,:));
        u_k(k)=mean(u_k_m(k,:));
    end
    
    clearvars a c ac
    
    cvx_begin quiet
    %cvx_begin
    variable xi_c
    variable P(Nt*(K+1),1) complex
    % variable a(K,K) complex
    % variable c(K)
    % variable ac(K,K) complex
    % variable cc(K)
    % PP=reshape(P,Nt,K+1);
    % PK=PP(:,2:K+1);
    % Pc=PP(:,1);
    for k=1:K
        for i=1:K
            a(k,i)=P((Nt*i+1):Nt*(i+1),1)'*psi_k(:,:,k)*P((Nt*i+1):Nt*(i+1),1);
            ac(k,i)=P((Nt*i+1):Nt*(i+1),1)'*psi_ck(:,:,k)*P((Nt*i+1):Nt*(i+1),1);
        end
        c(k)=real(f_k(:,k)'*P((Nt*k+1):Nt*(k+1),1));
    end
    
    
    xx=xi_c+sum(sum(a))+sum(N0*t_k+u_k-v_k)+sum(-2.*c);
    maximize(-xx)
    subject to
    for k=1:K
        P(1:Nt,1)'*psi_ck(:,:,k)*P(1:Nt,1)+ sum(ac(k,:))+N0*t_ck(k)-2*real(f_ck(:,k)'*P(1:Nt,1))+u_ck(k)-v_ck(k)<=xi_c;
    end
     P(1:Nt,1)==zeros(Nt,1);
    P'*P<=Pt;
    cvx_end
    
    %%%
    PP=reshape(P,Nt,K+1);
    pc=PP(:,1);
    pk=PP(:,2:K+1);
    
    n=n+1;
    A(n)=xx;
%     P'*P
%     Pt
    if A(n)>A(n-1)&&n>3
        fprintf('NOT monotonically decreasing %d in MMSE-RS \n',A(n)-A(n-1))
        A
    end
end
A;
end