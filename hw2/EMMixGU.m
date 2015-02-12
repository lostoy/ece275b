function [alpha1,mu,sigma,theta,Ps]=EMMixGU(X,theta,maxIter,opt)

N=size(X,1);
P=zeros(N,2);
Ps=zeros(1,maxIter);

[~,mus,Sigmas]=EMMixG(X,2,1000,0);
Sigma=min(Sigmas);
mu=mus(Sigma==Sigmas);

alpha1=0.5;

L0=-Inf;

theta0=theta;
tao=4;
beta=0.2;

for iter=1:maxIter
    
    P(:,1)=1/theta*(X<=theta&X>=0);
    P(P(:,1)==0 & X>=0)=1/theta0/100;
    %P(:,1)=1/(2*theta)*exp(-X/2/theta);
    P(:,2)=normpdf(X,mu,sqrt(Sigma));
    
    P=P*diag([alpha1,1-alpha1]);
    
    %current Log likelihood
    L=sum(log(sum(P,2)));
    
    %stopping criteira
    if ((L-L0)<0.001*abs(L0) && iter>30)
        break;
    end;
    L0=L;
    
    %true log likelihood
    Ps(iter)=L;
    
    %E-step
    Phat=diag(1./sum(P,2))*P;
    
    %M-step
    
    if (opt)
         theta1=Phat(:,1)'*X*2/sum(Phat(:,1));
    
         ind=P(:,1)>1/(10*theta1);
         theta2=max(X(ind).*exp(-(X(ind)-theta1).^2/tao^2));
         
         if isempty(theta2)
             theta2=theta1;
         end
         %[theta1,theta2]
         theta=beta*theta1+(1-beta)*theta2;
          
            
    else
         theta=max(X(Phat(:,1)>0));
    end
    
    
    
    mu=Phat(:,2)'*X/sum(Phat(:,2));
    Sigma=sum((X-repmat(mu,N,1)).^2.*Phat(:,2))/sum(Phat(:,2));
    alpha1=sum(Phat(:,1))/N;
    
end
Ps=Ps(1:iter-1);
sigma=sqrt(Sigma);
end

