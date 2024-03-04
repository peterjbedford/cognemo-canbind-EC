function [T_pr,score0] = cognemo_SVMtest(mdl,X_te)
%% Preamble
%{
%}
%% predict values and scores

[T_pr,score] = predict(mdl,X_te);
score0 = score(:,1);

end