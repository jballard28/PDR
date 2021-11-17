function [f,df]=rank_one(theta,~,~)

f=@(t)(t'*theta).^2;

% == 2*theta'*(t*t')
df=@(t)2*t'.*(t'*theta);
end

