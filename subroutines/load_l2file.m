function [snp,data,column_names,noSNPs] = load_l2file(filepath,suffix)
%LOAD_L2FILE loads a file of LD scores

expected_SNPs = 1e6;
chr=cell(22,1);
snp=chr;pos=chr;cm=chr;data=chr;


counter = 0;
noSNPs = 0;
for i=1:22
    l2file=fopen([filepath,num2str(i),suffix]);
    line=strsplit(fgets(l2file));
    [~,column_ind,columns_found] = intersect(line,...
        {'CHR','SNP','BP','CM'},'stable');
    if ~any(columns_found==2)
        error('SNP column not found')
    end
    
    fields = {'%d ', '%*c%*c%d ', '%d ', '%f '};
    matching_str = horzcat(fields{columns_found}, repmat('%f ',1,length(line)-length(column_ind)-1), '\n');
    
    column_names = line(length(column_ind):end-1);
    
    l2data=textscan(l2file,matching_str,'TreatAsEmpty',{'NA'});
    
    if any(columns_found==1)
        chr{i}=l2data{column_ind(columns_found==1)};
    end
    
    snp{i}=l2data{column_ind(columns_found==2)};
    
    if any(columns_found==3)
        pos{i}=l2data{column_ind(columns_found==3)};
    end
    if any(columns_found==4)
        cm{i}=l2data{column_ind(columns_found==4)};
    end
    data{i}=[l2data{length(column_ind)+1:length(line)-1}];
    fclose(l2file);
    
    if nargout > 3
        noSNPs = noSNPs + load([filepath,num2str(i),'.l2.M']);
    end
end

snp=vertcat(snp{:});
data=vertcat(data{:});
