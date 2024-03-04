function mdl = cognemo_SVMtrain(X,T,options)
%% Preamble
%{
%}
%% unpack options

% number of trees in forest
kernel = 'linear'; % DEFAULT
if isfield(options,'kernel')
    kernel = options.kernel;
end

%% train model

mdl = fitcsvm(X,T,...
                'KernelFunction',kernel);
%{
mdl = fitclinear(X,T,...
                'Learner','logistic',...
                'Regularization','ridge');
%}
end