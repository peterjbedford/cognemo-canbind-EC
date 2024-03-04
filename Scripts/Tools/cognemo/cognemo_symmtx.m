function Xsym = cognemo_symmtx(X,method)
%% Preamble
%{
Turns vectorized matrix into a symmetrical vectorized matrix
%}
%%
N_o = size(X,1);
N_c = size(X,2); N_r = sqrt(N_c);

Xsym = zeros(N_o,N_c);

for i = 1:N_o
    Xi = X(i,:);
    Xmtx = reshape(Xi,[N_r,N_r]);
    Xsymmtx = Xmtx + Xmtx';
    if method == 1
        % for when input is a triu or tril
        Xsym(i,:) = reshape(Xsymmtx,[1,N_c]);
    elseif method == 2
        % for when input is an asymmetrical matrix--gives average across triu,
        % tril
        Xsym(i,:) = reshape(Xsymmtx./2,[1,N_c]);
    else
        Xsym(i,:) = Xi;
    end
end
end