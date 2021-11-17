function [f, df]=uniform_corr(x,~,~)

f=@(t)((1-x)*sum(t.^2,1)+x*sum(t,1).^2)';
df=@(t)(sum(t,1).^2-sum(t.^2,1))';
end