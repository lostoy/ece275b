% %%10.a
% N=1000;
% alpha=0.55;
% mu=3/2;
% sigma=1/4;
% theta=1;
% 
% %histogramming and pdf
% X=randMixGU(alpha,mu,sigma,theta,N);
% [ph,x]=hist(X,20);
% ph=ph/(x(2)-x(1))/N;
% pdf=MixGUPdf(x,alpha,mu,sigma,theta);
% figure
% plot(x,ph,':*',x,pdf,'-s','LineWidth',2);
% title('histogram and pdf of Mixture Guassian and Uniform','FontSize',15,'FontWeight','Bold');
% legend('histogram','pdf')
% set(gca,'FontSize',15,'FontWeight','Bold');
% %TOSAVE
% %% MixG
% for K=[1 2 3 10 50]
% [alphas,mus,Sigmas,Ps]=EMMixG(X,K,10000);
% 
% XG=randMixG(alphas,mus,Sigmas,N);
% [phG,xg]=hist(XG,20);
% phG=phG/(xg(2)-xg(1))/N;
% pdf=MixGUPdf(xg,alpha,mu,sigma,theta);
% figure
% plot(xg,phG,':*',xg,pdf,'-s','LineWidth',2);
% title(['histogram of learned GMM(K=' num2str(K) ') and pdf of MGU'],'FontSize',15,'FontWeight','Bold');
% legend('histogram of GMM','pdf of MGU')
% set(gca,'FontSize',15,'FontWeight','Bold');
% %TOSAVE
% 
% figure
% plot(1:length(Ps),Ps,'LineWidth',2);
% title(['Loglikelihood K=',num2str(K)],'FontSize',15,'FontWeight','Bold');
% set(gca,'FontSize',15,'FontWeight','Bold');
% end
%% 

%% 10.a
for mu=[2 3/2 1 1/2]
    sol10a(mu);
end

%% 10.b
clear;clc;
N=2000;
alpha1=0.55;

sigma=1/4;
theta=1;

for mu=[2 3/2 1 1/2 ]

%sampling
X=randMixGU(alpha1,mu,sigma,theta,N);

%%% UE learning
% 
% [alphas,mus,sigmas,thetas,Ps]=EMMixGUE(X,max(X)+1,1000);
% thetas
% alphas
% figure;
% plot(1:length(Ps),Ps,'LineWidth',2);
% hold on;
% figure;
% x=-abs(mu)-1:0.01:abs(mu)+1;
% pdf=MixGUPdf(x,alpha1,mu,sigma,theta);
% pdfs=MixGUPdf(x,alphas,mus,sigmas,thetas);
% plot(x,pdfs,x,pdf);

%%% learning
[alphas,mus,sigmas,thetas,Ps]=EMMixGU(X,max(X)+1,1000,1);
thetas
alphas
figure;
plot(1:length(Ps),Ps,'LineWidth',2);
title(['Likelood \mu=' num2str(mu)],'FontSize',15,'FontWeight','Bold');
set(gca,'FontSize',15,'FontWeight','Bold');
        
xlabel('iteration');
ylabel('likelihood');
saveas(gca, ['./eps/10b/L' num2str(mu) '.eps'] ,'epsc');

figure;
x=-abs(mu)-1:0.1:max(theta+1,abs(mu)+1);
pdf=MixGUPdf(x,alpha1,mu,sigma,theta);
pdfs=MixGUPdf(x,alphas,mus,sigmas,thetas);
plot(x,pdfs,'--*',x,pdf,'-s','LineWidth',2);
title(['pdf \mu=' num2str(mu)],'FontSize',15,'FontWeight','Bold');
set(gca,'FontSize',15,'FontWeight','Bold');
legend('EM','groundtruth')

saveas(gca, ['./eps/10b/pdf' num2str(mu) '.eps'] ,'epsc');

end


%%
% the=6;
% x=[-1:0.1:the];
% 
% plot(x,1/0.5924/the*exp(-x/0.5924/the)*0.5+1/0.5924/the*exp(-(the-x)/0.5924/the)*0.5);
% ylim([0,1])

%%
% [alphas,mus,sigmas,thetas,Ps]=EMMixGEE(X,max(X)+1,1000);
% thetas
% alphas
% figure;
% plot(1:length(Ps),Ps,'LineWidth',2);
% hold on;
% 
% 
% x=-abs(mu)-1:0.01:abs(mu)+1;
% pdf=MixGUPdf(x,alpha1,mu,sigma,theta);
% [y1,y2]=EEpdf(x,thetas);
% figure;
% plot(x,alphas*0.5*(y1+y2)+(1-alphas)*normpdf(x,mus,sigmas),x,pdf);
