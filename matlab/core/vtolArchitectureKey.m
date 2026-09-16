function key = vtolArchitectureKey(architecture)
%VTOLARCHITECTUREKEY Normalize supported architecture names.
%
%   key = vtolArchitectureKey(architecture) returns either 'twin' or
%   'tail'.  Keeping this conversion in one place prevents a silent change
%   in actuator ordering between the model, allocator, and Simulink wrapper.

if nargin < 1 || isempty(architecture)
    error('vtolArchitectureKey:MissingArchitecture', ...
        'Specify ''twin''/''A'' or ''tail''/''B''.');
end

name = lower(strtrim(char(architecture)));
switch name
    case {'twin', 'twin-only', 'a', 'architecture a'}
        key = 'twin';
    case {'tail', 'tail-assisted', 'tri', 'tri-prop', 'b', 'architecture b'}
        key = 'tail';
    otherwise
        error('vtolArchitectureKey:UnknownArchitecture', ...
            'Unknown architecture ''%s''.', char(architecture));
end
end
