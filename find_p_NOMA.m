function [p_max]=find_p_NOMA(K,Pt,A,B,C,D,X,h)

if norm(h(:,1))>=norm(h(:,2))
    user1=1;
    user2=2;
else
    user1=2;
    user2=1;
end
p=[0;0;0;0;0;0];
[U,W,Z] = svds(X);%V=U*W*Z'
n=length(W);
obj_max=-100;
for m=1:50000
    r=sqrt(1/2)*(randn(n,1)+1i*randn(n,1));%sqrt(var/2)*(randn(1,N)+1i*randn(1,N))
    p1=U*sqrt(W)*r;
    p1(3:4)=[0;0];
    %     p2=[0;0;0;0;0;0];
    %     p2(user1:user1+1)=p1(1:2);
    %     p2(5:6)=p1(3:4);
    p=sqrt(Pt)*p1/norm(p1);
    p(3:4)=[0;0];
    p=sqrt(Pt)*p/norm(p);
    if (norm(h(:,1)'*p(1:2))>=norm(h(:,1)'*p(5:6)))&&(norm(h(:,2)'*p(1:2))>=norm(h(:,2)'*p(5:6)))
        
        [GMI,GMI_c,GMI_p]=cal_GMI(K,A,B,C,D,p);
        obj(m)=GMI;
        
        if obj(m)>=obj_max
            obj_max=obj(m);
            lc=min(GMI_c);
            add=GMI_p;
            p_max=p;
        end
    end
end

[obj_max,lc,add]*log(2);
%result==GMI*log(2)
end