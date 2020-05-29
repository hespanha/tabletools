function [scaledValues,scale]=scaleValues(rawValues,method,range)
% [scaledValues,scale]=scaleValues(rawValues,method,range)
% Automatically scales a vector
%
% [scaledValues,scale]=scaleValues(rawValues,scale)
% Scales a vector for a given power and coefficients
%
% Inputs:
%   rawValues - matrix of values to be scaled
%   method    - method used to automatically pre-process values 
%                'linear'            - no pre-processing
%                'perceptualDisks'   - assumes that raw values will be represented 
%                                      as areas of circles and returns values for
%                                      the radii that are "perceptually correct"
%                'perceptualSquares' - assumes that raw values will be represented 
%                                      as areas of squares and returns values for
%                                      the sides of the squares that are
%                                      "perceptually correct"
%                'mathematicalArea'  - assumes that raw values will be represented 
%                                      as areas of shapes and returns values for
%                                      scaling the shapes (corresponds to taking sqrt).
%   range []      - no scaling (a=1, b=0)
%     or
%   range [1 x 2] - desired range of values [minimum,maximum]
%     or
%   range [1 x 1] - desired maximum value (keeping 0 at 0)
%
%   scale  - structure with the scaling parameters:
%                scaledValues = (scale.a) * rawValues^(scale.p) + (scale.b)
%
% Outputs:
%   scaledValues - matrix with scaled values
%   scale        - structure with the scaling parameters
%                  (can be used as input to subsequent calls to scaleValues)
%
% Copyright (C) 2013-2014  Stacy & Joao Hespanha

if nargin==3 
    switch method
      case 'linear'
        p=1;
      case 'mathematicalArea'
        p=.5;
      case 'perceptualSquares'
        p=.54;
      case 'perceptualDisks'
        p=.57;
      otherwise
        error('unkown scaling method ''%s''\n',method)
    end
    scaledValues=rawValues.^p;

    if isempty(range)
        a=1;
        b=0;
    else
        mx=max(scaledValues(:));mn=min(scaledValues(:));
        
        if length(range)<2
            range=[0,range];
            mn=0;
        end
        
        if mx>mn
            a=(range(2)-range(1))/(mx-mn);
            b=range(1)-(range(2)-range(1))/(mx-mn)*mn;
        else
            fprintf('   ''range'' field is invalid when all values are equal (min=%g, max=%g), setting single value to range(2)\n',mn,mx);
            a=0;
            b=range(2);
        end
    end
else
    p=method.p;
    a=method.a;
    b=method.b;
    scaledValues=rawValues.^p;
end

scaledValues=a*scaledValues+b;
scale.a=a;
scale.b=b;
scale.p=p;

