function [f,df]=trait_specific(x,~,whichTrait)
f=@(t)x*t(whichTrait,:)'.^2;
df=@(t)t(whichTrait,:)'.^2;
end

