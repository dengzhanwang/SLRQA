function str = char(q)
% CHAR Create character array (string).
% (Quaternion overloading of standard Matlab function.)

% Note: the Matlab char function converts arrays of numeric values into
% character strings. This is not what this function does, but the Matlab
% guidance on user-defined classes suggests writing a char function and
% a disp/display function. This advice has been followed.

% Copyright © 2005, 2008, 2010 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(1, 1), nargoutchk(0, 1) 

% There are three cases to be handled. The argument is one of: empty, a
% pure quaternion, a full quaternion.

if isempty(q)
    str = '[] quaternion'; % This case must be handled first because an
    return;                % empty quaternion is not scalar, and the
end                        % next check would fail.

if ~isscalar(q)
    error('char cannot handle a vector or a matrix.')
end

% TODO The format here should be dependent on the current setting of the
% format (using the format command). This can be obtained using f = get(0,
% 'format') but unfortunately it returns a string like 'Short' rather than
% a C-style format code. There doesn't seem to be a built-in way to map one
% to the other, so we would have to do it here, with a switch statement.

f = '%.4g'; % Control over the format of each numeric value.

if isempty(q.w)
   % The scalar part is empty, so we begin with the x component of q.
   
   str = [nullminus(q.x) comp2str(q.x) ' * I'];
else
    % There is a scalar part, so we start with that, and then add in the x
    % component of q.
   str = [nullminus(q.w)     comp2str(q.w) ' ' ...
          plusminus(q.x) ' ' comp2str(q.x) ' * I'];
end

str = [str ' ' plusminus(q.y) ' ' comp2str(q.y) ' * J'];
str = [str ' ' plusminus(q.z) ' ' comp2str(q.z) ' * K'];

    function S = comp2str(X)
        % Create a string representation of one component of the quaternion,
        % which may be numeric (real or complex) or a symbolic, or logical.
 
        if isnumeric(X)
            if isreal(X)
                S = num2str(abs(X), f);
            else
                S = ['(' num2str(X, f) ')']; % Complex numeric.
            end
        elseif islogical(X)
            % Logical values cannot be complex (try it, any way you try to
            % make a complex with logical elements, they get converted to
            % double). Hence we can display logicals using the same code as
            % for numerics. We do it separately in case in the future we
            % decide to display logical values as T/F or true/false.
            S = num2str(abs(X), f);
        elseif isa(X, 'sym')
            % We form a string here with surrounding parentheses in some
            % cases and without in others. If the string contains an
            % operator like + ^ or / we put the parentheses in place. If
            % the string is just a variable name, we don't need to,
            % although we don't test for that explicitly at present. This
            % may need review.
            S = char(X);
            if ~isempty(regexp(S, '[\^/+-]', 'once'))
                S = ['(' S ')'];
            end
        else
            error(['Cannot convert quaternion with ', class(q.x), ...
                   'elements to a character representation'])
        end
    end
end

function S = plusminus(X)
% Extracts the sign of X and returns the character '-', '+'. The sign
% function doesn't exist for non-numeric values (e.g. logical) which is why
% we check whether X is numeric.

if isnumeric(X) && sign(X) == -1
    S = '-';
else
    S = '+';
end
end

function S = nullminus(X)
% Returns a space or minus according to the sign of X.
if isnumeric(X) && sign(X) == -1
    S = '-';
else
    S = ' ';
end

end

% $Id: char.m 1113 2021-02-01 18:41:09Z sangwine $
