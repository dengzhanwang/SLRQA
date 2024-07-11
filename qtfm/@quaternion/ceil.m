function a = ceil(q)
% CEIL   Round towards plus infinity.
% (Quaternion overloading of standard Matlab function.)

% Copyright © 2006 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(1, 1), nargoutchk(0, 1) 

a = overload(mfilename, q);

end

% $Id: ceil.m 1113 2021-02-01 18:41:09Z sangwine $
