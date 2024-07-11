function test_inv
% Test code for the quaternion inv function.

% Copyright © 2011 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

disp('Testing inverse ...')

T = 1e-10;

% Test scalar inverse on a few real and complex quaternions.

for i=1:10
    q = randq;
    compare(q * inv(q), onesq, T, 'quaternion/inv failed test 1') %#ok<*MINV>
    b = complex(q, randq) .* randn;
    compare(b * inv(b), onesq, T, 'quaternion/inv failed test 2')
end

A = randq(10);
compare(A * inv(A), eyeq(10), T, 'quaternion/inv failed test 3')

B = complex(randq(10), randq(10)) .* randn(10);
compare(B * inv(B), eyeq(10), T, 'quaternion/inv failed test 3')

disp('Passed');

% $Id: test_inv.m 1113 2021-02-01 18:41:09Z sangwine $
