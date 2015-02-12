function X=randMixG(alphas,mus,Sigmas,N)
D=size(mus(1,:),2);
K=length(alphas);

X=randn(N,D);

Salphas=zeros(K,1);
Salphas(1)=alphas(1);

for k=2:K
    Salphas(k)=Salphas(k-1)+alphas(k);
end
for i=1:N
    k=firstSmalller(rand,Salphas);
    X(i,:)=mvnrnd(mus(k,:),Sigmas(:,:,k));
end
end

function k=firstSmalller(x,a)
for k=1:length(a)
    if(x<=a(k))
        break;
    end
end

end