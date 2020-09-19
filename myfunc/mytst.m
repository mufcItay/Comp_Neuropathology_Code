function [data] = mytst(cfg)
%This functions simulates data for a single TST experiment.

%% set params
% in this code - the probability to choose MB or MF strategy depends on
% the realability of the task. we start with 0.5 Pmb (set in the sim.m
% file) and it will be modulated in the code according to the outcomes of
% the actions.
Pmb = cfg.Pmb; % instead of w in shira's  code. the probability to choose MB 
alpha1      = cfg.alpha1;  %Learning Rate for first stage
alpha2      = cfg.alpha2;  %Learning rate for second stage
eta_SPE     = cfg.eta_SPE;     %Learning rate for the transition matrix
eta_RPE     = cfg.eta_RPE;     %Learning rate for the transition matrix
Xmf         = 0.5;                % The initial realibility of MF strategy
Xmb         = 0.5;                % The initial realibility of MB strategy
A_alpha     = cfg.A_alpha;
B_alpha     = cfg.B_alpha;
A_beta      = cfg.A_beta;
B_beta      = cfg.B_beta;
beta        = cfg.beta;    %temperature parameters (the degree of randomness  in the action–selection). 
lambda      = cfg.lambda;  %eligibility trace
rewards     = cfg.rewards;
maxRPE      = max(rewards) - min(rewards); % we set the max RPE according to possible max and min rewards 

%% set vars
Ntrials     = cfg.Ntrials; %number of trails
Nstates     = cfg.Nstates; %number of states
Nstages     = cfg.Nstages; %number of stages
actions     = cfg.actions; % number of possible actions
Pmb_array = zeros(Ntrials,1); % the probability to choose MB strategie 
Pmb_array(1) = Pmb; 
Xmf_array = zeros(Ntrials + 1,1); % the reliability of MF
Xmf_array(1) = Xmf; 
Xmb_array = zeros(Ntrials + 1,1); % the reliability of MB
Xmb_array(1) = Xmb;
reward      = zeros(Ntrials,1); %a vector containing all the reward 
Qmf         = zeros(actions,Nstates,Nstages); %a matrix containing the Q values
Qnet        = zeros(2,1); %a matrix containing the Q values
pA          = zeros(Nstages,Ntrials); %predected probabilty for each action in a given state 
expected    = [.8,.5;.5,.5]; %a matrix containing the objective probabilty for a reward for each action in stage 2 (first row if we choose a1 at stage 1, first coloumn if we choose a1 at stage 2) ;a5,a6))
tran_prob_true = cfg.trans_prob_true;
tran_prob   = ones(actions, Nstates); % the subject will learn and update the transition matrix during the experiment
uniform_probability = 1/Nstates; % the transition matrix begins with the same values for all states
tran_prob = tran_prob*uniform_probability; % the probability is 1/number of atates

data=struct();

%% Run experiment
for t = 1:Ntrials 
    %first stage choice
    Qmb(1)  = tran_prob(1,1)* max(Qmf(:,1,2))+ tran_prob(1,2)*max(Qmf(:,2,2)); % action 1 value
    Qmb(2)  = tran_prob(2,1) * max(Qmf(:,1,2))+ tran_prob(2,2)*max(Qmf(:,2,2)); % action 2 value
    
    Qnet(1) = Pmb_array(t)*Qmb(1)+(1-Pmb_array(t))*Qmf(1,1,1);
    Qnet(2) = Pmb_array(t)*Qmb(2)+(1-Pmb_array(t))*Qmf(2,1,1);
    
    pA(1,t) = exp(beta.*Qnet(1))./...
             (exp(beta.*Qnet(1))+ exp(beta.*Qnet(2))); %the probabilty to choose one of the two actions in state1 (stage1). consists of the Qvalue of a1 & a2. 
    a1      = randsample([1,2],1,true,[pA(1,t) 1-pA(1,t)]);      %select 1 with probability pA and 2 with(1-pA)
    
    %transition into second stage
    state   = randsample([1,2],1,true,[tran_prob_true(a1,1) 1-tran_prob_true(a1,1)]);
    tran    =(a1==state)*1; %1 for common and 0 for uncommon

    % state prediction error calculation according to Lee et al 2013
    % supplamentry page 19.
    spe = 1 - tran_prob(a1,state);
    % according to the same paper, the subject updates the transition
    % matrix in this parttern:
    tran_prob(a1,state) = tran_prob(a1,state) + eta_SPE*spe; %
    % the alternative state has to sum the probablity to 1. because the
    % states are 1 or 2, 3 minus the satae will give us the other state:
    % 3-1=2, 3-2=1:
    tran_prob(a1,3-state) = 1 - tran_prob(a1,state); % 
    
    %second stage choice
    pA(2,t) = exp(beta.*Qmf(1,state,2))./...
             (exp(beta.*Qmf(1,state,2))+ exp(beta.*Qmf(2,state,2))); %the probabilty to choose one of the two actions in the cuurent state (stage2). consists of the Qvalue of the two actions of that state. 
    a2      = randsample([1,2],1,true,[pA(2,t) 1-pA(2,t)]);    %select 1 with probability pA and 2 with(1-pA)
    
    %outcome
    reward(t)=randsample(rewards,1,true,[1-expected(state,a2) expected(state,a2)]); %did we get a reward? retrieval of the relevant prob. out of the expected matrix
    
    %update Qvals
    PE1 = (Qmf(a2,state,2)-Qmf(a1,1,1)); %updating the predection eror of the choosen action at the first stage
    PE2 = (reward(t)-Qmf(a2,state,2));   %updating the predection eror of the choosen action at the choosen state in the second stage
    
    %save PE ans statePE history (history(trial) = sign(PE))
    if PE1 == 0, RPEsigned = PE1; else, if PE1 > 0, RPEsigned = 2; else, RPEsigned = 1; end, end  
    if spe == 0, SPEsigned = spe; else, if spe > 0, SPEsigned = 2; else, SPEsigned = 1; end, end  
    
    % D is the arrray that contains all PE signs in each trial = the history of
    % the prediction errors. (it can be: 0= no PE, 1 = negative PV or 2=
    % positive PE). we need this in order to calculate the realibailty of
    % each strategy. state-PE refers to MB and PE1 refers to MF (which is
    % reward-PE in the paper)
    D_RPE(t) = RPEsigned;
    D_SPE(t) = SPEsigned;
    
    Qmf(a1,1,1)    =Qmf(a1,1,1)+lambda*alpha1*PE1+alpha1*PE2; %updating the Q value of the action at the first stage
    Qmf(a2,state,2)=Qmf(a2,state,2)+alpha2*PE2;               % updating the Q value of the choosen action at the choosen state in the second stage
    
    % set params for reliabilty updating, and call respected function
    % set Xmf params
    paramsmf.Xmf = Xmf_array(t);
    paramsmf.D = D_RPE;
    paramsmf.PE1 = PE1;     
    paramsmf.maxRPE = maxRPE;     
    paramsmf.eta_RPE = eta_RPE;
    %update realbility values:
    Xmf_array(t+1) = cfg.reliabilityFunc(paramsmf);
    % set Xmb params
    paramsmb.Xmb = Xmb_array(t);
    paramsmb.D = D_SPE; 
    Xmb_array(t+1) = cfg.reliabilityFunc(paramsmb);

    transition_alpha = A_alpha / (1 + exp(B_alpha*Xmf_array(t+1))); 
    transition_beta = A_beta / (1 + exp(B_beta*Xmb_array(t+1)));
    
    % finaly, we can calculate the probability to choose MB and MF
    % strategies (we acctually calculate only the prob to choose MB)
    Pmb_array(t+1) = Pmb_array(t) + transition_alpha * (1 - Pmb_array(t)) - transition_beta * (Pmb_array(t));
    
    %keep data per trial
    data.t(t,:)     =t;
    data.a1(t,:)    =a1;
    data.tran(t,:)  =tran;
    data.state(t,:) =state;
    data.a2(t,:)    =a2;
    data.reward(t,:)=reward(t);
    data.Pmb (t,:) = Pmb_array (t);
    
end
end

