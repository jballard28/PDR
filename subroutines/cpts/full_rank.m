function [f,df] = full_rank(theta,~,~)

% convert theta into upper-tringular matrix
theta=triuind(theta);

% == diag(t'*theta'*theta*t)
f=@(t)sum((theta*t).^2,1)';

% d/dX t'X'Xt == 2Xtt'
df=@(t)cell2mat(arrayfun(@(col){triuind(2*theta*t(:,col)*t(:,col)')'},(1:size(t,2))'));

end

