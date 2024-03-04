function pdata = cognemo_prepC(C,data)
%% Preamble
%{
%}
%% Prepare data for analysis

Xin = C.rdata.X; T = logical(data.T);

[~,N_c] = size(Xin);

if inputname(1) == "FC"
    % get only triu off-diag
    [Xin_tu,~,~,tu_f_ind,~,~] = cognemo_splitmtx(Xin); tu_ind = find(tu_f_ind);
    % eliminate all-zero connections
    [Xin_dz,~,nz_f_ind,~,~] = cognemo_dezero(Xin_tu);
else
    % eliminate all-zero connections
    [Xin_dz,~,nz_f_ind,~,~] = cognemo_dezero(Xin);
end

f_ind = nz_f_ind;
% data to be used in analysis
X = Xin_dz; N_v = size(X,2);
X0 = X(~T,:); X1 = X(T,:);

% 'full' vectors (i.e. all connections have an entry)
X0_f = zeros(length(find(~T)),N_c); X0_f(:,f_ind) = X0;
X1_f = zeros(length(find(T)),N_c);  X1_f(:,f_ind) = X1;

%% Means, stds

MX0 = mean(X0,1); MX1 = mean(X1,1);
SX0 = std(X0,1);  SX1 = std(X1,1);

% 'full' vectors (i.e. all connections have an entry)
MX0_f = zeros(1,N_c); MX0_f(f_ind) = MX0;
MX1_f = zeros(1,N_c); MX1_f(f_ind) = MX1;
SX0_f = zeros(1,N_c); SX0_f(f_ind) = SX0;
SX1_f = zeros(1,N_c); SX1_f(f_ind) = SX1;

% 'full' vectors, symmetrical for plotting
MX0_f_sym = cognemo_symmtx(MX0_f,1);
MX1_f_sym = cognemo_symmtx(MX1_f,1);
SX0_f_sym = cognemo_symmtx(SX0_f,1);
SX1_f_sym = cognemo_symmtx(SX1_f,1);

%% Package output

pdata.X = X; pdata.N_v = N_v; pdata.N_c = N_c;
pdata.f_ind = f_ind;
pdata.X0 = X0; pdata.X0_f = X0_f;
pdata.MX0 = MX0; pdata.MX0_f = MX0_f; pdata.MX0_f_sym = MX0_f_sym;
pdata.SX0 = SX0; pdata.SX0_f = SX0_f; pdata.SX0_f_sym = SX0_f_sym;
pdata.X1 = X1; pdata.X1_f = X1_f;
pdata.MX1 = MX1; pdata.MX1_f = MX1_f; pdata.MX1_f_sym = MX1_f_sym;
pdata.SX1 = SX1; pdata.SX1_f = SX1_f; pdata.SX1_f_sym = SX1_f_sym;

end
