% Compute the improvement ratio
function ir = IR(alpha_true,alpha_hat,alpha1,alpha2,whichTrait)
    if nargin < 5
        whichTrait = 1:size(alpha_true,2);
    end

    alpha_true = alpha_true(:,whichTrait);
    alpha_hat = alpha_hat(:,whichTrait);
    alpha1 = alpha1(:,whichTrait);
    alpha2 = alpha2(:,whichTrait);

    ir = (SSE(alpha_true,alpha_hat) - SSE(alpha_true, alpha1))/...
        (SSE(alpha_true, alpha_hat) - SSE(alpha_true, alpha2));

end