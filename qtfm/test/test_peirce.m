function test_peirce
% Test code for the Peirce decomposition.

% Copyright © 2019 Stephen J. Sangwine and Nicolas Le Bihan
% See the file : Copyright.m for further details.

disp('Testing Peirce decomposition ...');

T = 1e-10;

q = randq(5) .* randn(5);

[chi, d, a, b] = peirce(q);

compare(q, a + b, T, 'quaternion/peirce failed test 1.');

compare(sin(q), sin(chi) .* conj(d) + sin(conj(chi)) .* d, ...
                  T, 'quaternion/peirce failed test 2.');

disp('Passed');

end

% $Id: test_peirce.m 1113 2021-02-01 18:41:09Z sangwine $
