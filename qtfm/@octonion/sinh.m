function Y = sinh(X)
% SINH   Hyperbolic Sine.
% (Octonion overloading of standard Matlab function.)

% Copyright © 2006 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(1, 1), nargoutchk(0, 1)

if isreal(X)
    
    % X is a real octonion, and we compute the hyperbolic sine of an
    % isomorphic complex number using the standard Matlab sinh function,
    % then construct an octonion with the same axis as the original
    % octonion.
    
    Y = isooctonion(sinh(isocomplex(X)), X);
else
    
    % X is a complex octonion, and therefore we cannot use the method
    % above for real octonions, because it is not possible to construct
    % an isomorphic complex number. We use instead a fundamental formula
    % for the hyperbolic sine.
    
    Y = (exp(X) - exp(-X)) ./ 2;
end

end

% $Id: sinh.m 1113 2021-02-01 18:41:09Z sangwine $

% Created automatically from the quaternion
% function of the same name on 17-Feb-2020.