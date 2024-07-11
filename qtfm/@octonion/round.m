function a = round(q)
% ROUND Round towards nearest integer.
% (Octonion overloading of standard Matlab function.)

% Copyright © 2006 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(1, 1), nargoutchk(0, 1) 

a = overload(mfilename, q);

end

% $Id: round.m 1113 2021-02-01 18:41:09Z sangwine $


% Created automatically from the quaternion
% function of the same name on 17-Feb-2020.