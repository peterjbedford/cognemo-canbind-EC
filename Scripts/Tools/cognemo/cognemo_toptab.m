function cognemo_toptab(x,xlabel,method,xset,inputname)
%% Preamble
%{
Converts input vector x (t-stat of difference, feat importance, etc.) which
is a 1xN_c (N_c:=number of connections) vector  representing connected
regions (i.e. N_c = N_r*N_r) to a table which displays the entries
rank-ordered by the values in x, labeled according to rlabel.

% rank #    region-from     region-to     x value

%}
%% 
N_top = xset.rankt.N_top;
rlabel = xset.rankt.rlabel;

N_c = length(x); N_r = sqrt(N_c);
rankN = string(1:N_top); % rank order

varnames = {'rank','region_from','region_to',xlabel};

%% Prep

[~,xtop_ind] = cognemo_topk(x,N_top,method);
xtop = x(xtop_ind);

[row,col]=ind2sub([N_r,N_r],xtop_ind); % get matrix index values

% region labels
r_from = rlabel(col); r_to = rlabel(row); 

T = table(rankN',r_from',r_to',xtop','VariableNames',varnames);

%% Prep ex

ex.x = x;
ex.x(xset.area_ex.c_ind) = 0;

[~,ex.xtop_ind] = cognemo_topk(ex.x,N_top,method);
ex.xtop = ex.x(ex.xtop_ind);

[ex.row,ex.col]=ind2sub([N_r,N_r],ex.xtop_ind); % get matrix index values

% region labels
ex.r_from = rlabel(ex.col); ex.r_to = rlabel(ex.row); 

ex.T = table(rankN',ex.r_from',ex.r_to',ex.xtop','VariableNames',varnames);

%% Write tables

rtdir = 'Ranking tables';
if not(isfolder('Results')); mkdir('Results'); end; cd 'Results'
if not(isfolder(rtdir)); mkdir(rtdir); end; cd(rtdir)
if not(isfolder(inputname)); mkdir(inputname); end; cd(inputname)

% check if already exists; delete and re-save
savename = [inputname '_' xlabel '_top_' char(string(N_top)) '.xlsx'];
if isfile(savename); delete(savename); end
writetable(T,savename);

%% Write tables ex

% check if already exists; delete and re-save
ex.savename = [inputname '_' xlabel '_top_' char(string(N_top)) '_ex'];
for i = 1:length(xset.area_ex.label)
    ex.savename = [ex.savename char(string(xset.area_ex.label(i)))];
    if i < length(xset.area_ex.label)
        ex.savename = [ex.savename '+'];
    end
end
ex.savename = [ex.savename '.xlsx'];

if isfile(ex.savename); delete(ex.savename); end
writetable(ex.T,ex.savename);

cd ../../..

end
