function [DCM,output,options] = cognemo_getEC_rDCM(Y,alloptions)
%% Preamble
%{
This function performs rDCM on a prepped BOLD data structure Y.
---------------------------------------------------------------------------
INPUTS
---------------------------------------------------------------------------
Y:=             the fMRI data, in a form required for rDCM analysis
alloptions:=     structure which includes:
- options:=     options for rDCM analysis
- est_method:=  chooses either regular or sparsity rDCM
- Ast:=         (only for method 1) a prior on the connectivity matrix,
                probably based on structural data
---------------------------------------------------------------------------
OUTPUTS
---------------------------------------------------------------------------
DCM:=           the DCM structure
output:=        the parameters as calculated through rDCM
options:=       options for rDCM analysis
%}

%%

options_in = alloptions.options;
options_in.y_dt = Y.dt;
est_method = alloptions.est_method;

% specify DCM
[DCM] = tapas_rdcm_model_specification(Y,[],[]);
data_type = 'r'; % empirical data

if isfield(alloptions,'Ast')
    DCM.a = alloptions.Ast;
    DCM.a(logical(eye(size(DCM.a)))) = 1;
end

% Estimate (with timer)
currentTimer = tic;
[output, options] = tapas_rdcm_estimate(DCM,data_type,options_in,est_method);
toc(currentTimer)
end