function [r]=R_square(yhat,y)

avg=mean(y);
r=1-sum((y-yhat).^2)/sum((y-avg).^2);


end