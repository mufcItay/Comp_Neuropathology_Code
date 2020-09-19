function X = bayesian_realability(params)
% this function is based on the supplemtary of the paper:
% "Neural Computations Underlying Arbitration between Model-Based and
% Model-free Learning"
% Sang Wan Lee, Shinsuke Shimojo, and John P. O’Doherty
% The calcultion of the realibilty will be according to: 
% "Bayesian Reliability estimation of MB and MF strategy" part of the paper

% input:
% the usual Params.
% we will use Params.D - an array that contains the history of PE type. 
% each place of the array contains a number that can be 1, 2 or, 0  (0= no 
% PE, 1 = negative PE or 2= positive PE) that indicates the type of PE that
% was in that trial (the index is the trial number) 

% output:
% the realability of the strategy. meaning Xmf or Xmb, depands which D we
% use when we call the function. in the mytst.m file we modulate params.D
% before we call this function - if we want to calculate Xmf we set params.D
% to be the history of PE calculated using MF strategy and vice versa.

% according to the paper we need to calculate the cardinality of D:
% from math-wiki: The cardinality of a set is a measure of a set's size, 
% meaning  = the number of elements in the set.
% in the case of this code, it refers to the number of the occurence of the
% events that leads to PEi (it can be: 0= no PE, 1 = negative PV or 2=
% positive PE). D = {PE(1),PE(2),..,PE(T)}. and the cardinality of D is the
% sum of the number of the occurence of the events that leads to PEi, for
% every i. simply: the cardinality of D = the size D.

D_cardinality = numel(params.D); % the size of D

% the Xmf or Xmb of the bayesian method are calculated as the ratio btwn 
% the expectation to some kind of PE (zero, neg or pos) and the varience
% of it.
% in the next part we calaulate the expextation and var of each PE type 
% according to D.

% first the number of the occurence of the events that lead to PEi (in the 
% paper it is singed as #PEi)
zero_PE_count = sum(params.D(:) == 0); 
neg_PE_count = sum(params.D(:) == 1);
pos_PE_count = sum(params.D(:) == 2);
% expectation calculation (according to the paper): 
zero_PE_exp = (1 + zero_PE_count) / (3 + D_cardinality); 
neg_PE_exp = (1 + neg_PE_count) / (3 + D_cardinality);
pos_PE_exp = (1 + pos_PE_count) / (3 + D_cardinality);
% variance calculation (according to the paper):
zero_PE_var = ((1 + zero_PE_count)*(2 + neg_PE_count + pos_PE_count))  / ((3 + D_cardinality)^2 * (4 + D_cardinality)); 
neg_PE_var = ((1 + neg_PE_count)*(2 + zero_PE_count + pos_PE_count))  / ((3 + D_cardinality)^2 * (4 + D_cardinality)); 
pos_PE_var = ((1 + pos_PE_count)*(2 + neg_PE_count + zero_PE_count))  / ((3 + D_cardinality)^2 * (4 + D_cardinality)); 
% now we calculate the X (realability) of each PE type:
X0 = zero_PE_exp/zero_PE_var;
X1 = neg_PE_exp/neg_PE_var;
X2 = pos_PE_exp/pos_PE_var;

% We now quantify the reliability of the learning strategy as the
% following (again according to the paper):
X = X0/(X0 + X1 +X2);

% that's it :)

end