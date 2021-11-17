function [SNPs,Z,phase] = import_sumstat_files(sumstatsPaths)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
nt=length(sumstatsPaths);
[SNPs,Z,phase] = import_sumstats(sumstatsPaths{1});
for ii=2:nt
    [newSNPs,newZ,newPhase] = import_sumstats(sumstatsPaths{ii});
    [SNPs,i1,i2]=intersect(SNPs,newSNPs,'stable');
    Z=[Z(i1,:),newZ(i2)]; %implement better
    phase=[phase(i1,:),newPhase(i2)];
end

    function [SNPs1,Z1,phase] = import_sumstats(filepath)
        %UNTITLED2 Summary of this function goes here
        %   Detailed explanation goes here
        file=fopen(filepath);
        fgets(file);
        data=textscan(file,'rs%d %s %s %*f %f %f %*f %*d %*d');
        if ~feof(file)
            error('Failed to load summary stats')
        end
        SNPs1=data{1};
        phase=(-1).^(cellfun(@(x,y)mod(sum(x(1)<y(1)),2)==0,data{2},data{3}));
        
        Z1=data{5}.*phase;
    end

end