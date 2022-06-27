function [p_max]=find_p(K,Pt,A,B,C,D,X,CC)

%[U,W,Z] = svds(X);%V=U*W*Z'

%U(:,1);%eigenvector for largest eigen value
%p_e=U(:,1)*sqrt(trace(X));

    [U,W,Z] = svds(X);%V=U*W*Z'
    n=length(W);
    obj_max=-100;
    for m=1:10000
        r=sqrt(1/2)*(randn(n,1)+1i*randn(n,1));%sqrt(var/2)*(randn(1,N)+1i*randn(1,N))
        p1=U*sqrt(W)*r;
        p=sqrt(Pt)*p1/norm(p1);
        obj(m)=cal_GMI_min(K,A,B,C,D,p,CC);
        if obj(m)>=obj_max
            obj_max=obj(m);
            p_max=p;
        end
    end

end