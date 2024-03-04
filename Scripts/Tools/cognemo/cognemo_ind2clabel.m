function clabel = cognemo_ind2clabel(ind,names,options)
%% Preamble
%{
Applies combined roi labels to indices ind to give connection labels.
---------------------------------------------------------------------------
INPUTS
---------------------------------------------------------------------------
ind:=               indices in linear form (1,N_v)
names:=             string array of names, raw from Y.name
options.exstr:=     a string which is to be stripped from all of the names
options.length:=    "short"=only allowed characters for tables etc;
                    else=no change
options.dir:=       1 if EC. Applies right arrow
---------------------------------------------------------------------------
OUTPUTS
---------------------------------------------------------------------------
clabel:=            the names of the connections indicated by the indices
                    ind
%}

%% Unpack options
middle = "\leftrightarrow";
if options.dir
    middle = "\rightarrow";
end
if isfield(options,'exstr')
    takeout = options.exstr;
end
if isfield(options,'length') && options.length == "short"
    takeout = [takeout " " "," "-"];
end
%%
N_r = length(names);
n = length(ind);

% convert indices from vector to indices in matrix
[row,col]=ind2sub([N_r,N_r],ind);

roirow = string(ones(1,n)); roicol = string(ones(1,n));


%%
for i_n = 1:n
    if isfield(options,'exstr') || (isfield(options,'length') && options.length == "short")
        roirow(i_n) = erase(string(names(:,row(i_n))),takeout);
        roicol(i_n) = erase(string(names(:,col(i_n))),takeout);
    else
        roirow(i_n) = string(names(:,row(i_n)));
        roicol(i_n) = string(names(:,col(i_n)));
    end
    if isfield(options,'areacolors')
        from_color = options.areacolors(row(i_n),:)./255;
        to_color   = options.areacolors(col(i_n),:)./255;
        from_colorstr = char(string(from_color(1)) + "," + ...
                             string(from_color(2)) + "," + ...
                             string(from_color(3)));
        to_colorstr   = char(string(to_color(1)) + "," + ...
                             string(to_color(2)) + "," + ...
                             string(to_color(3)));
        roirow(i_n) = "{\color[rgb]{"+from_colorstr+"}\bullet}"+roirow(i_n);
        roicol(i_n) = "{\color[rgb]{"+to_colorstr+"}\bullet}"+roicol(i_n);
    end
    if isfield(options,'length') && options.length == "short"
        roirowb4bkt = split(roirow(i_n),"(");
        roicolb4bkt = split(roicol(i_n),"(");
        roirow(i_n) = roirowb4bkt(1);
        roicol(i_n) = roicolb4bkt(1);
    end
end
clabel = roicol + middle + roirow;
clabel = erase(clabel," ");

end