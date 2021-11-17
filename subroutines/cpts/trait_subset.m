function [f,df] = trait_subset(theta,~,whichTrait)

% convert theta into upper-tringular matrix
theta=triuind(theta);

% == diag(t'*theta'*theta*t)
f=@(t)sum((theta*t(whichTrait,:)).^2,1)';


df=@(t)cell2mat(arrayfun(@(col){triuind(2*theta*t(whichTrait,col)*t(whichTrait,col)')'},(1:size(t,2))'));

end

