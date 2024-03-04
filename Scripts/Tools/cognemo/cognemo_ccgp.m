function cognemo_ccgp(C,xset)
%% Preamble
%{
Does the second set of connectogram plots:
- feature importance (averaged across direction for EC)
%}
%% PREP & SORT
N_c = C.pdata.N_c;
method = C.xset.symmethod;

CIMP = zeros([1,N_c]);
CIMP(C.pdata.f_ind) = mean(C.class.ccout.IMP,1);
[CIMPt,~] = cognemo_topk(CIMP,xset.cgp.N_top,2);
Cdata = ...
    cognemo_cg_prep(cognemo_symmtx(CIMPt,method),xset.cgp.ind);

clear CIMPt CIMP

%% Prep

% custom colormap 
N_r = sqrt(N_c);
Ccolormap = repmat(C.xset.pcolor,[N_r,1]);

%% Navigate to folder

cdir = 'Classification connectogram plots';
if not(isfolder('Results')); mkdir('Results'); end; cd 'Results'
if not(isfolder(cdir)); mkdir(cdir); end; cd(cdir)
if not(isfolder(C.inputname)); mkdir(C.inputname); end; cd(C.inputname)

%% Plot

Cplot = figure('Visible','off');
Cplot.Position = [50 50 xset.cgp.sx xset.cgp.sy];
circularGraphfi(Cdata,'Label',cellstr(xset.cgp.label),...
                     'ColorMap',Ccolormap);
T = cognemo_ringlabel;

%% Save and exit folder

saveas(Cplot,[C.inputname '_imp.png'])
cd ..\..\..

end