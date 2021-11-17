function [output] = triuind(input,symmetric)
%Indices of the upper-triangular entries of a square matrix of size N

if isscalar(input) && round(input)==input && input>0
    ind=arrayfun(@(n){(1:n)+input*(n-1)},1:input);
    output=[ind{:}]';
elseif isvector(input)
    N=floor(sqrt(2*length(input)));
    if N*(N+1)/2 ~= length(input)
        error('Input vector should have length N(N+1)/2 for some integer N')
    end
    output=zeros(N);
    output(triuind(N))=input;
    if exist('symmetric')
        if symmetric
            output = output + output' - diag(diag(output));
        end
    end
elseif size(input,1)==size(input,2)
    N=length(input);
    output=input(triuind(N));
else
    error('Input to triuind should either be a positive integer, a square matrix, or a vector')
end

end

