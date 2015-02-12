close all;clc;clear
alphas=[0.6 0.4];
mus=[0;10];
Sigmas=zeros(1,1,2);

allSigma=[1,25];
for Sigma=allSigma
    Sigmas(1,1,1)=Sigma;
    Sigmas(1,1,2)=Sigma;
    
    
    for N=[25 50 100 200 400];
        
        X=randMixG(alphas,mus,Sigmas,N);
        
        %histogram
        x=linspace(-10,20,20);
        [ph]=hist(X,x);
        ph=ph/(x(2)-x(1))/N;
        
        %pdf
        pdf=MixGPdf(x',alphas,mus,Sigmas);
        
        figure;
        plot(x,ph,':*',x,pdf,'-s','LineWidth',2);
        title(['histogram and pdf of Mixture Guassians \sigma^2=' num2str(Sigma) ' N=' num2str(N)],'FontSize',15,'FontWeight','Bold');
        
        set(gca,'FontSize',15,'FontWeight','Bold');
        
        %
        
        
        [alphams,mums,Sigmams,Ps]=EMMixG(X,2,1000,0);
        [alphams1,mums1,Sigmams1,Ps1]=EMMixG(X,2,1000,1);
        hold on
        
        plot(x,MixGPdf(x',alphams,mums,Sigmams),':hk','LineWidth',2);
        hold on;
        plot(x,MixGPdf(x',alphams1,mums1,Sigmams1),'--+c','LineWidth',2);
        title(['pdf of Mixture Guassians \sigma^2=' num2str(Sigma) ' N=' num2str(N)],'FontSize',15,'FontWeight','Bold');
        legend('histogram','pdf','EM','EM(same sigma)')
        set(gca,'FontSize',15,'FontWeight','Bold');
        
        if N==400
            
            mu_bar=mean(X);
            Sigma_bar=1/N*sum(X.^2)-mu_bar.^2;
            
            mu_bar_true=alphas(1)*mus(1)+alphas(2)*mus(2);
            Sigma_bar_true=alphas(1)*(mus(1)^2+Sigmas(1))+alphas(2)*(mus(2)^2+Sigmas(2))-mu_bar_true^2;
            hold on;
            plot(x,normpdf(x,mu_bar,sqrt(Sigma_bar)),'-^m','LineWidth',2);
            legend('histogram','pdf','EM','EM(same sigma)','one gaussian')
            disp(['mu--------Sigma']) 
            disp([mu_bar,Sigma_bar,mu_bar_true,Sigma_bar_true]);
        end
        
        saveas(gca, ['./eps/8/' num2str(N) '_' num2str(Sigma) '.eps'] ,'epsc');

        figure;
        plot(1:length(Ps),Ps,'LineWidth',2);
        title(['Likelood \sigma^2=' num2str(Sigma) ' N=' num2str(N)],'FontSize',15,'FontWeight','Bold');
        set(gca,'FontSize',15,'FontWeight','Bold');
        
        xlabel('iteration');
        ylabel('likelihood');
        
        saveas(gca, ['./eps/8/P' num2str(N) '_' num2str(Sigma) '.eps'] ,'epsc');

    end
end