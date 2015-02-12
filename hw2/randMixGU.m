function [x]=randMixGU(alpha,mu,sigma,theta,N)

x=zeros(N,1);
for i=1:N
    if rand<alpha
        x(i)=rand*theta;
    else
        x(i)=randn*sigma+mu;
        
    end
end
end