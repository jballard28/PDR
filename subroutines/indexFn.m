function [indexArray] = indexFn(linInds, noChoices)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here


indexArray = zeros(length(linInds), length(noChoices));
for i = 1:numel(noChoices)
    indexArray(:,i) = mod(linInds,noChoices(i));
    linInds = (linInds - indexArray(:,i)) / noChoices(i);

end

