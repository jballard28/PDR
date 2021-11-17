% Compute the sum of squared errors; this function calculates the sse for
% each column separately, then adds them afterwards
function sse = SSE(alpha1,alpha2)

    for i=1:size(alpha1,2)
        sse(i) = sum((alpha1(:,i)-alpha2(:,i)).^2);
    end
    
    sse = sum(sse);

end