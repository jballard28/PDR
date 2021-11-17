function [f, df]=null_cpt(~,~,~)
f=@(t)zeros(size(t,2),1);
df=@(t)zeros(size(t,2),0);
end