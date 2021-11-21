figure;
set(gcf,'Position',[249 492 1485 420]);

% set number of pleiotropic components
npleio = 2;

tnames = {'T2D','BMI','TG'};

totcov = est.cov;
newcpts = {est.cpts.cov};

subplotn = [1 2 4 5];

for ii=1:npleio
    
    cptnum = est.noCpts-npleio+ii;
    h_temp = newcpts{cptnum};
    D = diag(1./sqrt(diag(totcov)));
    h_norm = D*h_temp*D;
    
    subplot(1,npleio+1,ii);
    
    h1 = heatmap(tnames,tnames,h_norm);
    h1.FontSize = 16;
    h1.CellLabelFormat = '%.2f';
    colorbar off;
    colormap(bluewhitered(256))
    caxis([-1 1])
    h{ii} = h_norm;
    
end

sbp = subplot(1,npleio+1,npleio+1);
sbp.Position = sbp.Position .* [1.1 1 0.2 1];

matsum = sum(cat(3,h{:}),3);
t_sp = ones(size(matsum,1),1)-diag(matsum);
h1 = heatmap('_',tnames,t_sp);
h1.FontSize = 16;
h1.CellLabelFormat = '%.2f';
colormap(bluewhitered(256))
caxis([-1 1])
title('trait-specific components')