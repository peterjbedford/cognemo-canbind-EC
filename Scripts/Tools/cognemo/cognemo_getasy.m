function asydata = cognemo_getasy(C)
%% Preamble
%{
%}
%%
N_c = C.pdata.N_c;
N_r = sqrt(N_c);

% triu
u_ind  = logical(reshape(triu(ones(N_r),1),[1,N_c]));
X0u_f  = zeros(size(C.pdata.X0_f)); X0u_f(:,u_ind) = C.pdata.X0_f(:,u_ind);
X1u_f  = zeros(size(C.pdata.X1_f)); X1u_f(:,u_ind) = C.pdata.X1_f(:,u_ind);
MX0u_f = zeros(1,N_c); MX0u_f(u_ind) = C.pdata.MX0_f(u_ind);
MX1u_f = zeros(1,N_c); MX1u_f(u_ind) = C.pdata.MX1_f(u_ind);

% tril
l_ind  = logical(reshape(tril(ones(N_r),-1),[1,N_c]));
X0l_f  = zeros(size(C.pdata.X0_f)); X0l_f(:,l_ind) = C.pdata.X0_f(:,l_ind);
X1l_f  = zeros(size(C.pdata.X1_f)); X1l_f(:,l_ind) = C.pdata.X1_f(:,l_ind);
MX0l_f = zeros(1,N_c); MX0l_f(l_ind) = C.pdata.MX0_f(l_ind);
MX1l_f = zeros(1,N_c); MX1l_f(l_ind) = C.pdata.MX1_f(l_ind);


% align indices for direct comparison
X0l_f_sym = cognemo_symmtx(X0l_f,1); X0l_f_sym(:,l_ind) = 0; X0l_f = X0l_f_sym;
X1l_f_sym = cognemo_symmtx(X1l_f,1); X1l_f_sym(:,l_ind) = 0; X1l_f = X1l_f_sym;

% disregard connections which are zero for both triu and tril
X0_f = zeros(size(X0l_f,1)+size(X0u_f,1),N_c);
X1_f = zeros(size(X1l_f,1)+size(X1u_f,1),N_c);
indu = 1:size(X0u_f,1);
indl = (size(X0u_f,1)+1):(size(X0l_f,1)+size(X0u_f,1));
X0_f(indu,:) = X0u_f; X0_f(indl,:) = X0l_f;
X1_f(indu,:) = X1u_f; X1_f(indl,:) = X1l_f;
[X0,~,nz_f_ind0,~,~] = cognemo_dezero(X0_f);
[X1,~,nz_f_ind1,~,~] = cognemo_dezero(X1_f);
X0u = X0(indu,:); X0l = X0(indl,:);
X1u = X1(indu,:); X1l = X1(indl,:);

%% Prepare output

asydata.f_ind0 = nz_f_ind0; asydata.f_ind1 = nz_f_ind1;

asydata.u_ind = u_ind;
asydata.X0u = X0u; asydata.MX0u_f = MX0u_f;
asydata.X1u = X1u; asydata.MX1u_f = MX1u_f;

asydata.l_ind = l_ind;
asydata.X0l = X0l; asydata.MX0l_f = MX0l_f;
asydata.X1l = X1l; asydata.MX1l_f = MX1l_f;

end
