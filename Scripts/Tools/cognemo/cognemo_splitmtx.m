function [X_tu,X_tl,X_di,tu_ind,tl_ind,di_ind] = cognemo_splitmtx(X)
%% Preamble
%{
Has to be before de-zeroing
%}
%%
N_c = size(X,2); N_r = sqrt(N_c);

% construct indices, for convenience
tu_mtx = triu(ones(N_r),1);  tu_ind = logical(reshape(tu_mtx,[1,N_c]));
tl_mtx = tril(ones(N_r),-1); tl_ind = logical(reshape(tl_mtx,[1,N_c]));
di_mtx = eye(N_r);           di_ind = logical(reshape(di_mtx,[1,N_c]));

% divide data
X_tu = zeros(size(X)); X_tu(:,tu_ind) = X(:,tu_ind);
X_tl = zeros(size(X)); X_tl(:,tl_ind) = X(:,tl_ind);
X_di = zeros(size(X)); X_di(:,di_ind) = X(:,di_ind);

end
    