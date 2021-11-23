function plot_heatmap(est,tnames,npleio)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

totcov = est.cov;
newcpts = {est.cpts.cov};

% determining configuration of subplots
if npleio < 4
    subplotrow = 1;
    subplotcol = npleio+1;
    subplotn = 1:npleio+1;
elseif mod(npleio+1,3) == 2
    nrow = floor((npleio+1)/3);
    subplotn = [];
    for rr = 1:nrow-1
        subplotn = [subplotn (1:3)+ones(1,3)*3*(rr-1)];
    end
    subplotn = [subplotn (1:2)+ones(1,2)*3*(nrow-1) (1:3)+ones(1,3)*3*nrow];
    subplotrow = nrow+1;
    subplotcol = 3;
elseif mod(npleio+1,3) == 1
    subplotrow = npleio/3;
    subplotcol = 4;
    subplotn = [];
    for rr = 1:subplotrow-1
        subplotn = [subplotn (1:3)+ones(1,3)*4*(rr-1)];
    end
    subplotn = [subplotn (1:4)+ones(1,4)*4*(subplotrow-1)];
else
    subplotrow = (npleio+1)/3;
    subplotcol = 3;
    subplotn = 1:npleio+1;
end


for ii=1:npleio
    
    cptnum = est.noCpts-npleio+ii;
    h_temp = newcpts{cptnum};
    D = diag(1./sqrt(diag(totcov)));
    h_norm = D*h_temp*D;
    
    
    subplot(subplotrow,subplotcol,subplotn(ii));
    
    
    h1 = heatmap(tnames,tnames,h_norm);
    h1.FontSize = 16;
    h1.CellLabelFormat = '%.2f';
    colorbar off;
    colormap(bluewhitered(256))
    caxis([-1 1])
    h{ii} = h_norm;
    
end

sbp = subplot(subplotrow,subplotcol,subplotn(end));
sbp.Position = sbp.Position .* [1.1 1 0.2 1];

matsum = sum(cat(3,h{:}),3);
t_sp = ones(size(matsum,1),1)-diag(matsum);
h1 = heatmap('_',tnames,t_sp);
h1.FontSize = 16;
h1.CellLabelFormat = '%.2f';
colormap(bluewhitered(256))
caxis([-1 1])
title('trait-specific components')

end

