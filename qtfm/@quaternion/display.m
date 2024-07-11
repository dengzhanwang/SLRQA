function display(q) %#ok<DISPLAY>
% DISPLAY Display array.
% (Quaternion overloading of standard Matlab function.)

% Copyright © 2005, 2008 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

if isempty(inputname(1))
    name = []; % Changed behaviour 8 June 2020 to match what MATLAB now
               % does, which is to display the value alone in the case
               % where there is no named variable. Previously it would show
               % 'ans = ' before the value.
else
    name = inputname(1);
end

% The test for Matlab on the next line is a crude hack to make this code
% work under Octave, since Octave (3.2.4 at least) does not recognise the
% FormatSpacing parameter.
if ismatlab && isequal(get(0,'FormatSpacing'),'compact')
    if ~isempty(name)
        disp([name ' =']);
    end
    disp(q);
else
    if ~isempty(name)
        disp(' ');
        disp([name ' =']);
        disp(' ');
    end
    disp(q);
    disp(' ');
end
end

function TF = ismatlab
% Returns true if the code is running under Matlab.

S = ver('Matlab');
TF = ~isempty(S);
end

% $Id: display.m 1113 2021-02-01 18:41:09Z sangwine $
