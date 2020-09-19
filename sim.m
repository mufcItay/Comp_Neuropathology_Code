%examining the main effect& interaction while changing the W and the alpha 
close all

%set vars

cfg.trans_prob_true = [0.9 0.1; 0.1 0.9];
cfg.Ntrials = 1000;  %number of trails according to rndwlk
cfg.Nstates = 2;    %number of states
cfg.Nstages = 2;    %number of stages
cfg.actions = 2;    % number of possible actions
cfg.rewards = [0 1];% the possible rewards
Nsub        = 50;

res         =table();
for v1=0:0.2:1

disp([v1]);

for subj=1:Nsub

%prealocating
if subj==1; data2=table(); end 

%set params
cfg.alpha1  = v1;                    %Learning Rate first stage
cfg.alpha2  = rand(1);               %Learning Rate second stage

%% choose between Pearce-Hall and Bayesian reliability updating metohds (itay and adva)
default_Pmb = 0.5;
%cfg.reliabilityFunc = @PH_reliability;
cfg.reliabilityFunc = @bayesian_realability;

%% parameters that are taken from Lees supp paper (itay and adva)
% maybe there more vars that we implemented, but here is the part that gets
% complicated so we seperated it, please see page 23 in the paper to understand:

cfg.eta_SPE     = rand(1)/100;       %Learning rate for the transition matrix
cfg.eta_RPE     = rand(1)/100;       %Learning rate for the MF realibility

%% Model Based params (0.1,0.3)
% cfg.alphaX1 = 0.1; % transition rate from MF to MB of the first trial(1)in which Xmf = 1, is a fixed value, see page 23 in Lees's supp paper
% cfg.A_alpha = 0.3; % has to be bigger than: 2 * cfg.alphaX1; % maximum transition rate from MF to MB
%% Model Free params (0.05,0.1)
cfg.alphaX1 = 0.05; % transition rate from MF to MB of the first trial(1)in which Xmf = 0.1, is a fixed value, see page 23 in Lees's supp paper
cfg.A_alpha = 0.1; % has to be bigger than: 2 * cfg.alphaX1; % maximum transition rate from MF to MB
cfg.B_alpha = log(cfg.alphaX1^-1 * cfg.A_alpha -1); % the steepness rate to transform from MF to MB
% transition rate from MB to MF = beta: 
%% Model Based params (0.01,0.03)
% cfg.betaX1 = 0.01; % transition rate from MB to MF of the first trial(1) in which Xmf = 0.01, is a fixed value, see page 23 in Lees's supp paper
% cfg.A_beta = 0.03; % has to be bigger than: 2 * cfg.betaX1; % maximum transition rate from MB to MF
%% Model Free params (0.1,0.3)
cfg.betaX1 = 0.1; % transition rate from MB to MF of the first trial (1) in which Xmf = 0.01, is a fixed value, see page 23 in Lees's supp paper
cfg.A_beta = 0.3; % has to be bigger than: 2 * cfg.betaX1; % maximum transition rate from MB to MF
cfg.B_beta = log(cfg.betaX1^-1 * cfg.A_beta -1); % the steepness rate to transform from MB to MF
cfg.Pmb = default_Pmb;
%%
cfg.beta    = myrandrng(1,0.1,10);   %temperature parameters (the degree of randomness  in the action–selection). 
cfg.lambda  = rand(1);               %eligibility trace

%sim experiment & add relevant columns
[data]          = mytst(cfg);

data.reward_prev= mylag(data.reward,1)';
data.tran_prev  = mylag(data.tran,1)';
data.pStay      = (data.a1==mylag(data.a1,1)')*1;

%calculate for prob (CR, UR, CU, UU)
prob=[mean(data.pStay(data.reward_prev==1 &data.tran_prev==1));
      mean(data.pStay(data.reward_prev==1 &data.tran_prev==0));
      mean(data.pStay(data.reward_prev==0 &data.tran_prev==1));
      mean(data.pStay(data.reward_prev==0 &data.tran_prev==0))];

%aggrgate data across subjects
data2=[data2;table(cfg.alpha1,prob(1),prob(2),prob(3),prob(4),...
         'VariableNames',{'alpha1','CR', 'UR', 'CU', 'UU'})];
end
res=[res;data2];
end

bar(prob)
ylim([0,1])
res.maineff=(res.CR+res.UR)./2-(res.CU+res.UU)./2;
res.intereff=(res.CR-res.UR)-(res.CU-res.UU);
xticklabels({'R common' 'R rare' 'U common' 'U rare'});
