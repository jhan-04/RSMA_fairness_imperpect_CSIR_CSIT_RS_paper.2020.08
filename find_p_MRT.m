function [p_max]=find_p_MRT(M,Pt,A,B,C,D,Xc,pk)

[U,W,Z] = svds(Xc);%V=U*W*Z'
n=length(W);
obj_max=-100;
for m=1:10000
    r=sqrt(1/2)*(randn(n,1)+1i*randn(n,1));%sqrt(var/2)*(randn(1,N)+1i*randn(1,N))
    pc1=U*sqrt(W)*r;
    pc=sqrt(trace(Xc))*pc1/sqrt(pc1'*pc1);
    p=[pk;pc];
%     for k=1:M
%         GMI_p(k)=log2((p'*C(:,:,k)*p)/(p'*D(:,:,k)*p));
%         GMI_c(k)=log2((p'*A(:,:,k)*p)/(p'*B(:,:,k)*p));
%     end
%     
    obj(m)=cal_GMI(M,A,B,C,D,p);
    
    if obj(m)>=obj_max
        obj_max=obj(m);
%         lc=min(GMI_c);
%         add=sum(GMI_p);
        p_max=p;
    end
end

%[obj_max,lc,add]*log(2);
%result==GMI*log(2)
end