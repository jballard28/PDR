    function [f, df]=uniform_corr_slope(x,noTraits,whichTrait) 
        % x(1) corresponds to the correlation
        % x(2:end) corresponds to the diagonal elements other than
        %   whichTrait
        
        v=ones(noTraits,1);
        v([1:whichTrait-1,whichTrait+1:end])=x(2:end);
        
        temp = ones(noTraits)*x(1);
        temp = temp + eye(noTraits)*(1-x(1));
        
        S=diag(v)*temp*diag(v);
        
%         f=@(times)arrayfun(@(t)t'*S*t, times);
        f=@(times)splitapply(@(t)t'*S*t, times, 1:size(times,2))';
        
        h=1e-30;
        x = x+h;
        
        v=ones(noTraits,1);
        v([1:whichTrait-1,whichTrait+1:end])=x(2:end);
        temp = ones(noTraits)*x(1);
        temp = temp + eye(noTraits)*(1-x(1));
        S1=diag(v)*temp*diag(v);
        
        f1=@(times)splitapply(@(t)t'*S1*t, times, 1:size(times,2))';
        
        df=@(t)(f1(t)-f(t))/h;
        
    end
    
% function [S]=uniform_corr_slope(x,noTraits,whichTrait) 
%     % x(1) corresponds to the correlation
%     % x(2:end) corresponds to the diagonal elements other than
%     %   whichTrait
%     error('Not yet implemented')
%     v=ones(noTraits,1);
%     v([1:whichTrait-1,whichTrait+1:end])=x(2:end);
%     S=diag(v)*uniform_corr(x(1),noTraits)*diag(v);
% 
% end