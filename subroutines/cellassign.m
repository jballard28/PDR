function [cellArray] = cellassign(values,cells)
%CELLASSIGN takes a vector of length n and a cell array of vectors whose
%size sums to n, and it assigns the new values to the cells (preserving
%size of each cell)
indices=[0,cumsum(cellfun(@numel,cells))];
indices=arrayfun(@(i0,i1){i0+1:i1},indices(1:end-1),indices(2:end));
cellArray=cellfun(@(ii){values(ii)},indices);
end

