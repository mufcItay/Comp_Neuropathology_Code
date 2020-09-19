function X = PH_reliability(params)
% this function is based on the supplemtary of the paper:
% "Neural Computations Underlying Arbitration between Model-Based and
% Model-free Learning"
% Sang Wan Lee, Shinsuke Shimojo, and John P. O’Doherty.
% The calcultion of the realibilty will be according to: 
% "Pearce-Hall associability for reliability estimation" part of the paper

% input:
% the usual Params.

% output:
% the realability of the strategy. meaning Xmf or Xmb

% we can calculate using pearce hall (PH) calculation only the realability of MF
% strategy (and not Xmb). So we usuallu will use the baysian calculation instead. 
% if we are using PH we will return the same realability of MB that was
% before the estimation of this function. 

if(isfield(params, 'Xmb'))
    X = params.Xmb;
    return
end

%calculated according to the same paper Lee et al 2013 (Pearce-Hall associability for reliability estimation)
X = params.Xmf + params.eta_RPE *((1-abs(params.PE1)/params.maxRPE) -params.Xmf); 

% that's it :)

end