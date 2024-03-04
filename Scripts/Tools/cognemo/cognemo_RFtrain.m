function mdl = cognemo_RFtrain(X,T,options)
%% Preamble
%{
%}
%% unpack options

% number of trees in forest
N_tree = 100; % DEFAULT
if isfield(options,'N_tree')
    N_tree = options.N_tree;
end

%% train model

mdl = TreeBagger(N_tree,X,T,...
                'Method','classification',...
                'PredictorSelection','curvature',...
                'OOBPredictorImportance','on');
            
end