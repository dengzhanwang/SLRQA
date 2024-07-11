function tf = isnilpotent(q, tol)
% ISNILPOTENT  True where any element of q is a nilpotent to within the
% tolerance given (optionally) by the second parameter.

% Copyright © 2019 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(1, 2), nargoutchk(0, 1)

if nargin == 1
    tol = 8 .* eps; % The tolerance was not specified, supply a default.
end

% Reference:
%
% Stephen J. Sangwine and Daniel Alfsmann,
% Determination of the Biquaternion Divisors of Zero, Including the
% Idempotents and Nilpotents, Advances in Applied Clifford Algebras, 20,
% (2010), 401–410. DOI 10.1007/s00006-010-0202-3.

% Theorem 5 of the above paper gives the conditions for a biquaternion to
% be a nilpotent. The quaternion must be pure, and the norm must be zero.

if ispure(q)
    tf = abs(normo(q)) < tol; % The whole array is pure, test the norms.
else
    % Nilpotents must be pure, so we check that the scalar part is close to
    % zero, then check the norm of the vector part is also close to zero.
    
    tf = s(q) < tol & abs(normo(v(q))) < tol;
end

end

% $Id: isnilpotent.m 1113 2021-02-01 18:41:09Z sangwine $
