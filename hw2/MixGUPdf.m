function Ps=MixGUPdf(X,alpha1,mu,sigma,theta)

Ps=alpha1*normpdf(X,mu,sigma)+(1-alpha1)*1/theta*(X<theta&X>0);

end