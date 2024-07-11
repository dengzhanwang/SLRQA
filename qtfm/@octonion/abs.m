function a = abs(o)
% ABS Absolute value, or modulus, of an octonion.
% (Octonion overloading of standard Matlab function.)

% Copyright © 2011, 2015, 2021 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(1, 1), nargoutchk(0, 1)

if isnumeric(x(o.a))
    if isreal(o)
        
        % We use here a method based on Cayley-Dickson form and the Matlab
        % function hypot, which avoids overflow for large values.
        
        a = hypot(abs(o.a), abs(o.b));
    else
        
        % TODO The complexified octonion code is the same as the symbolic
        % and logical code. Merge them.
        
        % Here, for bioctonions, we cannot use hypot.
        
        a = sqrt(normq(o.a) + normq(o.b));
    end
    
elseif isa(x(o.a), 'sym') || islogical(x(o.a))

    a = sqrt(normq(o.a) + normq(o.b));

else
    error(['Unable to handle octonion with ' class(o.x) ' components'])
end

end

% $Id: abs.m 1125 2021-08-15 20:43:38Z sangwine $
