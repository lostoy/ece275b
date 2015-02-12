function [y1,y2]=EEpdf(x,theta)

b0=0.59245427;
y1=1/(theta*b0)*(exp(-x/(theta*b0)).*(x>=0));
y2=1/(theta*b0)*(exp(-(theta-x)/(theta*b0)).*(x<=theta &x>=0));
end