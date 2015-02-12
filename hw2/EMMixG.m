function [alphas,mus,Sigmas,Ps]=EMMixG(X,K,maxIter,option)


N=size(X,1);
D=size(X,2);

alphas=ones(1,K)/K;
[~,mus]=kmeans(X,K);
%mus=X(randperm(N,K),:);
Sigmas=repmat(eye(D),1,1,K);
Ps=zeros(maxIter,1);

L0=-Inf;

for iter=1:maxIter
    P=zeros(N,K);
    for k=1:K
        P(:,k)=mvnpdf(X,mus(k,:),Sigmas(:,:,k));
    end
    P=P*diag(alphas');
    %current Log likelihood
    L=sum(log(sum(P,2)));
    
    %stopping criteria
    if ((L-L0)<0.001*abs(L0) && iter>30)
        break;
    end
    L0=L;
    
    %True likelihood
    Ps(iter)=L;
    
    %E-step
    Phat=diag(1./sum(P,2))*P;
    %M-step
    sumtmp=0;
    for k=1:K
       mus(k,:)=Phat(:,k)'*X/sum(Phat(:,k));
       if option==0
           Sigmas(:,:,k)=(X-repmat(mus(k,:),N,1))'*diag(Phat(:,k))*(X-repmat(mus(k,:),N,1))/sum(Phat(:,k));
       else
           sumtmp=sumtmp+(X-repmat(mus(k,:),N,1))'*diag(Phat(:,k))*(X-repmat(mus(k,:),N,1));
       end
       alphas(k)=sum(Phat(:,k))/N;
    end
    if (option==1)
        for k=1:K
            Sigmas(:,:,k)=sumtmp/N;
        end;
    end
end
Ps=Ps(1:iter-1);
end