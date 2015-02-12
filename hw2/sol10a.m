function sol10a( mu )
N=1000;
alpha1=0.55;
sigma=1/4;
theta=1;

%histogramming and pdf
X=randMixGU(alpha1,mu,sigma,theta,N);
[ph,x]=hist(X,20);
ph=ph/(x(2)-x(1))/N;
pdf=MixGUPdf(x,alpha1,mu,sigma,theta);
figure
plot(x,ph,':*',x,pdf,'-s','LineWidth',2);
title(['histogram and pdf of MGU ( mu=' num2str(mu) ')'],'FontSize',15,'FontWeight','Bold');
legend('histogram','pdf')
set(gca,'FontSize',15,'FontWeight','Bold');
%TOSAVE
%% MixG
for K=[1 2 3 10 50]
[alphas,mus,Sigmas,Ps]=EMMixG(X,K,10000,0);

XG=randMixG(alphas,mus,Sigmas,N);
[phG,xg]=hist(XG,20);
phG=phG/(xg(2)-xg(1))/N;
pdf=MixGUPdf(xg,alpha1,mu,sigma,theta);
figure
plot(xg,phG,':*',xg,pdf,'-s','LineWidth',2);
title(['histogram of GMM(K=' num2str(K) ') and pdf of MGU ( mu=' num2str(mu) ')'],'FontSize',15,'FontWeight','Bold');
legend('histogram of GMM','pdf of MGU')
set(gca,'FontSize',15,'FontWeight','Bold');
%TOSAVE

figure
plot(1:length(Ps),Ps,'LineWidth',2);
title(['Loglikelihood K=',num2str(K),'mu=' num2str(mu)],'FontSize',15,'FontWeight','Bold');
set(gca,'FontSize',15,'FontWeight','Bold');
%TOSAVE
end

end

