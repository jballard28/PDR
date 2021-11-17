function [eigenvalues] = import_eig(filenames)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
eigenvalues = cell(1,1e4);
counter = 0;
for i = 1:numel(filenames)
    file = fopen(filenames{i});
    while ~feof(file)
        counter = counter + 1;
        line=split(fgets(file));
        eigenvalues{counter} = cellfun(@str2num,line(1:end-1));
    end
end
eigenvalues = eigenvalues(1:counter);
end

