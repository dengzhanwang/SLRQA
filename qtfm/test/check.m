function check(L, E)
% Test function to check a logical condition L, and output an error
% message from the string in the parameter E if false. L may be a vector,
% in which case, all its elements are required to be true.

% Copyright © 2005, 2008 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(2, 2), nargoutchk(0, 0)

if ~islogical(L)
    error('First parameter must be logical.');
end

if ~all(L)
    error(E);
end

% $Id: check.m 1113 2021-02-01 18:41:09Z sangwine $

