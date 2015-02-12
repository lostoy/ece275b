function Ps=MixGPdf(X,alphas,mus,Sigmas)

N=size(X,1);
K=length(alphas);
Ps=zeros(N,K);
for k=1:K
    Ps(:,k)=mvnpdf(X,mus(k,:),Sigmas(:,:,k));
end
Ps=Ps*alphas';
    
end