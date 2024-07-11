function x = sylvester(a, b, c)
% SYLVESTER  Solve quaternion sylvester equation ax + xb = c for x.
% (Quaternion overloading of standard Matlab function.)

% Copyright © 2016, 2020 Stephen J. Sangwine.
% See the file : Copyright.m for further details.

% Reference:
%
% Drahoslava Janovska and Gerhard Opfer,
% "Linear equations in quaternionic variables",
% Mitteilungen der Mathematischen Gesellschaft in Hamburg 27 (2008),
% 223–234. Theorem 2.3. Available: https://www.math.uni-hamburg.de/mathges/

narginchk(3, 3), nargoutchk(0, 1)

if ~isa(a, 'quaternion') || ~isa(b, 'quaternion') || ~isa(c, 'quaternion')
    error('All parameters must be quaternions.')
end

na = normq(a);
nb = normq(b);

if na <= nb
    % Use the second formula in Theorem 2.3 of the paper above.
    
    ib = b.^-1;
    
    x = (c + conj(a) .* c .* ib) .* (2 .* s(a) + b + na .* ib).^-1;
    
else
    % Use the first formula.
    
    ia = a.^-1;
    
    x = (2 .* s(b) + a + nb .* ia).^-1 .* (c + ia .* c .* conj(b));
end

end

% Note: An earlier version of this function released in QTFM 2.5 used
% adjoint matrices and the MATLAB sylvester function.

% $Id: sylvester.m 1113 2021-02-01 18:41:09Z sangwine $
