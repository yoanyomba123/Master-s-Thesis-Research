function Parameter = segparainitial(Parameter, imgTemplate)
dim1 = size(imgTemplate, 1); %#ok<*NASGU>
dim2 = size(imgTemplate, 2);
if ( size(Parameter.initialLSFBase,2) == 4 ) % which means you specify the initial LSF ...
  % by setting a rectangular Region of Interest.
    initialLSFCoor = Parameter.initialLSFBase;
    Parameter.initialLSF = ones([512 512]);
    Parameter.initialLSF(initialLSFCoor(2):initialLSFCoor(4),initialLSFCoor(1):initialLSFCoor(3)) = -1;
else
    if (imgTemplate == 0)
        error('You must specify the initial LSF by 4 Coors or other method.');
    else
        logicTemp = (imgTemplate == 0.021);
        initialLSF = ones(dim1, dim2);
        Parameter.initialLSF = logicTemp .* -2 + initialLSF;
    end
end