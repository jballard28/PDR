function [f,df]=fixed(~,~,Sigma)
f=@(t)sum(t.*(Sigma*t),1)';
df=@(t)zeros(size(t,2),0);
end

