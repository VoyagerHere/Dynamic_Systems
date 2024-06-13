function o = odeget(options,name,default,flag)
%ODEGET Get ODE OPTIONS parameters.
%   VAL = ODEGET(OPTIONS,'NAME') extracts the value of the named property
%   from integrator options structure OPTIONS, returning an empty matrix if
%   the property value is not specified in OPTIONS.  It is sufficient to
%   type only the leading characters that uniquely identify the
%   property.  Case is ignored for property names.  [] is a valid OPTIONS
%   argument.
%   
%   VAL = ODEGET(OPTIONS,'NAME',DEFAULT) extracts the named property as
%   above, but returns VAL = DEFAULT if the named property is not specified
%   in OPTIONS.  For example
%   
%     val = odeget(opts,'RelTol',1e-4);
%   
%   returns val = 1e-4 if the RelTol property is not specified in opts.
%   
%   See also ODESET, ODE45, ODE23, ODE113, ODE15S, ODE23S.

%   Mark W. Reichelt and Lawrence F. Shampine, 3/1/94
%   Copyright (c) 1984-98 by The MathWorks, Inc.
%   $Revision: 1.27 $  $Date: 2003/03/12 13:16:17 $
if (nargin==4) & isequal(flag,'fast')
    o=getknownfield(options,name,default);
    return
end

if nargin < 2
  error('Not enough input arguments.');
end
if nargin < 3
  default = [];
end

if ~isempty(options) & ~isa(options,'struct')
  error('First argument must be an options structure created with ODESET.');
end

if isempty(options)
  o = default;
  return;
end

Names = [
    'AbsTol      '
    'BDF         '
    'Events      '
    'InitialStep '
    'Jacobian    '
    'JacobianP   '
    'Hessians    '
    'HessiansP   '
    'Der3        '
    'Der4        '
    'Der5        '
    'JConstant   '
    'JPattern    '
    'Mass        '
    'MassConstant'                      % grandfathered
    'MassSingular'
    'MaxOrder    '
    'MaxStep     '
    'NormControl '
    'NonNegative '
    'OutputFcn   '
    'OutputSel   '
    'Refine      '
    'RelTol      '
    'Stats       '
    'Vectorized  '
    ];
[m,n] = size(Names);
names = lower(Names);

lowName = lower(name);
j = strmatch(lowName,names);
if isempty(j)               % if no matches
  error(sprintf(['Unrecognized property name ''%s''.  ' ...
                 'See ODESET for possibilities.'], name));
elseif length(j) > 1            % if more than one match
  % Check for any exact matches (in case any names are subsets of others)
  k = strmatch(lowName,names,'exact');
  if length(k) == 1
    j = k;
  else
    msg = sprintf('Ambiguous property name ''%s'' ', name);
    msg = [msg '(' deblank(Names(j(1),:))];
    for k = j(2:length(j))'
      msg = [msg ', ' deblank(Names(k,:))];
    end
    error('%s).', msg);
  end
end

if any(strcmp(fieldnames(options),deblank(Names(j,:))))
  eval(['o = options.' Names(j,:) ';']);
  if isempty(o)
    o = default;
  end
else
  o = default;
end
% --------------------------------------------------------------------------
function v = getknownfield(s, f, d)
%GETKNOWNFIELD  Get field f from struct s, or else yield default d.

if isfield(s,f)   % s could be empty.
  v = subsref(s, struct('type','.','subs',f));
  if isempty(v)
    v = d;
  end
else
  v = d;
end
