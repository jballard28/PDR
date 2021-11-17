    function [f,df]=uncorr(x,~,~)
        S=diag(x);
        
        f=@(t)sum(t.*(S*t),1)';
        df=@(t)t'.^2;
    end
