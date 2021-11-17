%% Fitting models

%% Setting up model to fit
no_traits = est_data.noTraits;
noBlocks = 100;

scalars=[0 5.^(0:4)];
h2Con='none';%'diagonal','full'

model_additive_cpts=1;
nonneg=1;

% setup true model
clear components;

for tt=1:no_traits
    S=zeros(no_traits);
    S(tt,tt)=1;
    components(tt) = COMPONENT('fixed',no_traits,scalars,S,...
        'h2ConOpt', h2Con);
end

for cc=1:4
    pleiotropic_cpt_type = 'full_rank';
    components(end+1) = COMPONENT(pleiotropic_cpt_type,no_traits,scalars,1,...
        'h2ConOpt', h2Con);
end

for cc=1:0
    pleiotropic_cpt_type = 'rank_one';
    components(end+1) = COMPONENT(pleiotropic_cpt_type,no_traits,scalars,1,...
        'h2ConOpt', h2Con);
end

h2ConVal=[];
est = MODEL(components, h2ConVal, model_additive_cpts);

%% ECF
samplingtimes=sqrt(1./scalars(2:end));
ecf=ECF(est_data,samplingtimes,'tol',.1,'noBlocks',noBlocks);

%% Fitting model

% ninit = 5000;
% nrot = 20;
% ninit = 8000;
% nrot = 20;
ninit = 10000;
nrot = 20;
% ninit = 1000;
% nrot = 20;

convergence_param=1e-6;

% gdsteps = 4000;
gdsteps = 6000;
% gdsteps = 1000;

est.initialize_fit(ecf,'niter',ninit,'method',"init_cov",'datacov',est_data.cov,...
    'gdsteps',gdsteps,'nrot',nrot,'convergence_param',convergence_param);

%% save variables

save(join(['realtraits_091021_4fullrank_',clustername]),'est_data','est','ecf');

beep