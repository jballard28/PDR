%% Fitting models

%% Setting up model to fit
no_traits = data.noTraits; % Number of traits, inferred from data object
noBlocks = 100; % Number of jackknife blocks for calculating empirical characteristic function (ECF)

scalars=[0 5.^(0:4)]; % Scalars for each mixture gaussian within each component
h2Con='none';%'diagonal','full'

model_additive_cpts=1;
nonneg=1;

% Set up model scaffold: number and type of components
clear components;

% Trait-specific components
% Generally, we have a trait-specific component for each trait
for tt=1:no_traits
    S=zeros(no_traits);
    S(tt,tt)=1;
    components(tt) = COMPONENT('fixed',no_traits,scalars,S,...
        'h2ConOpt', h2Con);
end

% Generic pleiotropic components (full-rank)
no_fr = 2; % How many generic pleiotropic cpts
for cc=1:no_fr
    pleiotropic_cpt_type = 'full_rank';
    components(end+1) = COMPONENT(pleiotropic_cpt_type,no_traits,scalars,1,...
        'h2ConOpt', h2Con);
end

% Factor-like pleiotropic components (rank-one)
no_r1 = 0; % How many factor-like pleiotropic components
for cc=1:no_r1
    pleiotropic_cpt_type = 'rank_one';
    components(end+1) = COMPONENT(pleiotropic_cpt_type,no_traits,scalars,1,...
        'h2ConOpt', h2Con);
end

h2ConVal=[];
est = MODEL(components, h2ConVal, model_additive_cpts);

%% ECF
samplingtimes=sqrt(1./scalars(2:end)); % Choosing sampling times based on the scalars
ecf=ECF(data,samplingtimes,'tol',.1,'noBlocks',noBlocks); % Caulculating the ECF

%% Fitting model

ninit = 100; % Number of initializations
nrot = 10; % Number of times to iteratively fit and update ECF

convergence_param=1e-6; % Parameter for indicating that convergence has been achieved

gdsteps = 5000; % Max number of gradient descent steps for parameter estimation step

est.initialize_fit(ecf,'niter',ninit,'method',"covmat",'datacov',data.cov,...
    'gdsteps',gdsteps,'nrot',nrot,'convergence_param',convergence_param,...
    'rotation_tolerance',1e-2);

%% save variables

save('fitmodel','data','est','ecf'); % Saves data object, model, and ECF to a file called fitmodel.mat

beep