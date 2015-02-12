function [alpha1,mu,sigma,theta,Ps]=EMMixGUE(X,theta,maxIter)

N=size(X,1);
%P=zeros(N,4);
P=zeros(N,3);

Ps=zeros(1,maxIter);

[~,mus,Sigmas]=EMMixG(X,2,1000);
Sigma=min(Sigmas);
mu=mus(Sigma==Sigmas);
%mu=X(randperm(N,1));
%Sigma=1;
alpha1=0.5;

L0=-Inf;


theta0=theta;

beta=0.98;
tao=4;
for iter=1:maxIter
     
     %[P(:,1),P(:,2),P(:,3)]=UEpdf(X,theta);
     [P(:,1),P(:,2),~]=UEpdf(X,theta);
     %P(:,4)=normpdf(X,mu,sqrt(Sigma));
      P(:,3)=normpdf(X,mu,sqrt(Sigma));  
     %P=P*diag([alpha1*0.33,alpha1*0.33,alpha1*0.33,1-alpha1]);
      P=P*diag([alpha1*0.99,alpha1*0.01,1-alpha1]);

        
    
    
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
    %mu=Phat(:,4)'*X/sum(Phat(:,4));
    %Sigma=sum((X-repmat(mu,N,1)).^2.*Phat(:,4))/sum(Phat(:,4));

    mu=Phat(:,3)'*X/sum(Phat(:,3));
    Sigma=sum((X-repmat(mu,N,1)).^2.*Phat(:,3))/sum(Phat(:,3));

    theta1=max(X(Phat(:,1)>0));
    theta2=Phat(:,2)'*(X-theta1)/(sum(Phat(:,1)+Phat(:,2)));
    
    theta=beta*theta1+(1-beta)*theta2;
    
    %theta1=(1-Phat(:,3))'*X*2/sum(1-Phat(:,3));
    
    %ind=(P(:,1)+P(:,2))>1/(10*theta1);
    
%    theta2=max(X(ind).*exp(-(X(ind)-theta1).^2/tao^2));
    
%    if isempty(theta2)
%         theta2=theta1;
%     end
%     theta=beta*theta1+(1-beta)*theta2;
    
    
    %M-alpha
%     P(:,1)=1/theta*(X<=theta&X>=0);
%     P(:,2)=normpdf(X,mu,sqrt(Sigma));
%     P(:,1:2)=P(:,1:2)*diag([alpha1,1-alpha1]);
%     
%     Phat=diag(1./sum(P(:,1:2),2))*P(:,1:2);
%     
%     alpha1=sum(Phat(:,1))/N;
    
    alpha1=sum(Phat(:,1)+Phat(:,2))/N;

end

Ps=Ps(1:iter-1);
sigma=sqrt(Sigma);

end
