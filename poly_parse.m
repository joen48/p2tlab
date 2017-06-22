function polysys = poly_parse(filename)
%POLY_PARSE Read set of polynomial equations from a file.
%   polysys = poly_parse(filename) reads set of polynomial equations from
%   filename and outputs the set as a cell.
%
% See also:
%   read_system

system = read_system(filename);
% number of equations
s = size(system, 1);
% variables
vars = {};
polysys = cell(s,2);
for i = 1:s
    S = sscanf(system{i,1}, '%s');
     A = regexp(S, '[0-9\.E\+\-]+(\*[A-Za-z0-9]+(\*\*[0-9]+)?)+(+|-)?', 'match');    % monomials
    p = size(A, 2);
    polysys{i,1} = zeros(p,1);
    polysys{i,2} = zeros(p,s);
    for l = 1:p 
        a = A{l};
        % coefficient
        c = str2double(regexp(a, '^[0-9\.E\+\-]+', 'match'));
        polysys{i,1}(l) = c;
        % variables
        x = regexp(a, '\*[A-Za-z0-9]+(\*\*[0-9]+)?', 'match');
        n = size(x, 2);
        for j = 1:n
            v = regexp(x{j}, '[A-Za-z0-9]+', 'match');
            v = v{:};
            if strcmpi(v, 'i')  % imaginary unit
                polysys{i,1}(l) = polysys{i,1}(l)*1i;
            else
                [LIA, LOCB] = ismember(v, vars);
                if ~LIA
                    vars{end+1} = v;
                    LOCB = size(vars, 2);
                end
                % exponent
                e = regexp(x{j}, '\*\*[0-9]+', 'match');
                if isempty(e)
                    e = 1;
                else
                    e = str2double(e{1}(3:end));
                end
                polysys{i,2}(l,LOCB) = e;
            end
        end
    end
end