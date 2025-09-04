function finder(varargin)
% opens a mac finder window in the path specified, or in the
% current directory if called without input
if nargin==0
    path = ['"', pwd, '"'];
else
    path = varargin{1};
end

eval(['! open ',path])
end