% Function to calculate piecewise linear function
function [betas, pred] = piecewise_linreg(x,y,thresholds)

% Making everything positive
signx = sign(x);
y = y.*signx;
x = x.*signx;

xs = zeros(length(x), length(thresholds));
for tt=1:length(thresholds)
    xs(:,tt) = min(x,thresholds(tt));
end

betas = xs\y;
pred = xs * betas;
pred = pred.*signx;
end