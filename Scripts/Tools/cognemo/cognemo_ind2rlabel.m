function rlabel = cognemo_ind2rlabel(ind,names,options)
%% Preamble
%{
extracts roi labels from Yfile, returns at given indices
---------------------------------------------------------------------------
INPUTS
---------------------------------------------------------------------------
ind:=               indices in linear form (1,N_v)
names:=             string array of names, raw from Y.name
options.exstr:=     a string which is to be stripped from all of the names
options.length:=    "short"=only allowed characters for tables etc;
                    else=no change
---------------------------------------------------------------------------
OUTPUTS
---------------------------------------------------------------------------
rlabel:=            the names of the connections indicated by the indices
                    ind
%}
%%
n = length(ind);

rlabel = string(ones(1,n));
takeout = options.exstr;
if isfield(options,'length') && options.length == "short"
    takeout = [takeout " " "," "-"];
end
for i_n = 1:n
    rlabel(i_n) = erase(names(:,ind(i_n)),takeout);
    if isfield(options,'length') && options.length == "short"
        roib4bkt = split(rlabel(i_n),"(");
        rlabel(i_n) = roib4bkt(1);
    end
end
rlabel = erase(rlabel," ");

end