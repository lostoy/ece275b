function [ y1,y2, y3] = UEpdf( x ,theta)

y1=1/theta*(x>=0 & x<=theta);
y2=1/(theta)*exp(-(x-theta)/(theta)).*(x>=theta);
y3=1/(0.5*theta)*exp(-x/(0.5*theta)).*(x>=0);
end

