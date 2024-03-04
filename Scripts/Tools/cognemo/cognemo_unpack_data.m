function fnames = cognemo_unpack_data(dataloc)
%% Preamble
%{
This function takes the content of the data folder and extracts the
filenames so as to begin the analysis.
---------------------------------------------------------------------------
INPUTS
---------------------------------------------------------------------------
dataloc:=   a folder name where the dataset files are contained

MATLAB
---------------------------------------------------------------------------
OUTPUTS
---------------------------------------------------------------------------
fnames:=    the list of filenames associated with the datasets
%}
%%
orig_dir = pwd;
cd(dataloc)
% Data file names

try
    fnames = string(ls('*.mat'));
    skip = startsWith(fnames,'.');
    idx_keep = find(~skip);
    fnames = fnames(idx_keep);
catch
    % this should be used for older versions of MATLAB
    fnames = ls;
    fnames = strsplit(fnames); N_fnames = length(fnames)-1;
    fnames = fnames(1:N_fnames);
    fnames = sort(fnames);
end

cd(orig_dir)
end
