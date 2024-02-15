 % [T,Y] = odeK4(@eqn,[время1,время2],[y0_1;y0_2],шаг,[42, 1],{d,[n1;n2],сигма,[F01;F02]},0,0);
[T,Y] = odeK4(@eqn,[10000,14000],[0; 0],0.001,[42, 1],{0.0005,[1;1],1,[0;0]},0,0);
    Y = mod(Y, 2*pi);
    
    plot(T, Y(:,1),'-r',T, Y(:,2),'-b')
    ylim([0 2*pi])
    yticks([0 pi 2*pi]); yticklabels(["0" "\pi" "2\pi"]);
title('d = 0.15', '\sigma = 1');
legend('n_1=1','n_2=1', 'Interpreter','tex');
xlabel("t");
ylabel('\theta', 'Interpreter','tex');


function [tout,xout,other] = odeK4(functionHandle,tspan,xinitial,timestep,extraparameters,parameters,next,quat)
%function [tout,xout] = odeK4(functionHandle,tspan,xinitial,timestep,extraparameters,next)
%%This function operates much like ode45 only it uses a fixed step
%RK4 integration
[N,flag] = size(xinitial);
if flag > 1
  disp('xinitial must be a column vector')
  tout = 0;
  xout = 0;
  return
end
tout = tspan(1):timestep:tspan(end);
integrationsteps = length(tout);
xout = zeros(integrationsteps,N);
x = xinitial;
if ~exist('extraparameters','var')
  extraparameters = 0;
  y = 0;
elseif isempty(extraparameters)
  extraparameters = 0;
  y = 0;
else
  y = extraparameters;
  extraparameters = 1;
  other = zeros(length(y),integrationsteps);
end
if ~exist('next','var')
  next = 0;
end
if strcmp(next,'off')
    next = 0;
end
threshold = 0;
if ~exist('quat','var')
  quat = 'none';
end
for ii = 1:length(tout)
  time = tout(ii);
  if next
    percent = (floor(100*(time-tout(1))/(tout(end)-tout(1))));
    if percent >= threshold
        
      disp(['Simulation ',num2str(percent),'% Complete'])
      
      threshold = threshold + next;
    end
  end
  %%Save States
  xout(ii,:) = x;
  %xdot(:,ii) = xdotRK4;
  %%Integrate
  
  if extraparameters
    if find(y == 42)
      [xdot1,outs] = feval(functionHandle, time, x, parameters{1}, parameters{2}, parameters{3}, parameters{4});
      other(:,ii) = outs;
      newF = outs;
      [xdot2,outs] = feval(functionHandle, time + (.5*timestep), x + (xdot1*.5*timestep), parameters{1}, parameters{2}, parameters{3}, parameters{4});
      [xdot3,outs] = feval(functionHandle, time + (.5*timestep), x + (xdot2*.5*timestep), parameters{1}, parameters{2}, parameters{3}, parameters{4});
      [xdot4,outs] = feval(functionHandle, time + timestep, x + (xdot3*timestep), parameters{1}, parameters{2}, parameters{3}, parameters{4});
      
      xdotRK4 = (1/6) * (xdot1 + (2*xdot2) + (2*xdot3) + xdot4);
      
      parameters{4} = newF;
      
      nextstate = x + (timestep * xdotRK4);
      
    else
      xdot1 = feval(functionHandle, time, x, y);
      xdot2 = feval(functionHandle, time + (.5*timestep), x + (xdot1*.5*timestep), y);
      xdot3 = feval(functionHandle, time + (.5*timestep), x + (xdot2*.5*timestep), y);
      xdot4 = feval(functionHandle, time + timestep, x + (xdot3*timestep), y);
      xdotRK4 = (1/6) * (xdot1 + (2*xdot2) + (2*xdot3) + xdot4);
      nextstate = x + (timestep * xdotRK4);
    end
  else
    xdot1 = feval(functionHandle, time, x);
    xdot2 = feval(functionHandle, time + (.5*timestep), x + (xdot1*.5*timestep));
    xdot3 = feval(functionHandle, time + (.5*timestep), x + (xdot2*.5*timestep));
    xdot4 = feval(functionHandle, time + timestep, x + (xdot3*timestep));
    xdotRK4 = (1/6) * (xdot1 + (2*xdot2) + (2*xdot3) + xdot4);
    nextstate = x + (timestep * xdotRK4);
  end
  if strcmp(quat,'Quat')
    %%Normalize the quaternions
    q0 = nextstate(4);
    q1 = nextstate(5);
    q2 = nextstate(6);
    q3 = nextstate(7);
    normq = sqrt(q0^2+q1^2+q2^2+q3^2);
    nexstate(4) = nextstate(4)/normq;
    nexstate(5) = nextstate(5)/normq;
    nexstate(6) = nextstate(6)/normq;
    nexstate(7) = nextstate(7)/normq;
  end
  x = nextstate;
end
end

function [dy_dt, F] = eqn(~, y, d, no, sigma, F)
  
  yy = mod(y, 2*pi);
  
  g = [1.001; 1.002];
  f = g - cos(y ./ no);
  exch = d * [F(2); F(1)];
  dy_dt = f - exch;
    
  [F1, F2] = chech_condition(yy, sigma);
  
  F = [F1; F2];
  
end

% Set F value on current phi value
function [F1, F2] = chech_condition(y, sigma)
  if ((y(1) > pi/2 - sigma) && (y(1) < pi/2 + sigma)) 
    F1 = 0;
  else
    F1 = 1;
  end

  if ((y(2) > pi/2 - sigma) && (y(2) < pi/2 + sigma)) 
    F2 = 0;
  else
    F2 = 1;
  end
end

function [DIFF, diff, ratio] = ms_phase_sync(T, Y)
    [A1, A2, err] = find_spikes(Y);
    if (err > 0) % death
        DIFF = NaN;
        diff = NaN;
        ratio = -err;
        return;
    end
    RATIO = zeros(length(A1(:, 1)) - 2, 1);
    for i = 2:length(A1) - 1
        RATIO(i - 1, 1) = sum(A2 < A1(i + 1, 1)) - sum(A2 < A1(i, 1));
    end
    ratio = mode(RATIO);
    NEAR = []; % find near spikes
    for i = 1:length(A1)
        [~, index] = max(A2(A2 <= A1(i)));
        if ~isempty(index)
            NEAR(end + 1) = A2(index);
        end
    end
    
    AA1 = id(A1, T);
    
    NEAR_ID = id(NEAR, T);
    len = length(NEAR_ID(:, 1));
    DIFF = zeros(len, 1);
    AA1 = AA1((end - len + 1):end, 1);
    DIFF(:, 1) = AA1 - NEAR_ID;
    diff = DIFF(end, 1);
end

function [A1, A2, err] = find_spikes(Y)
    error = 0.05; % Tolerance for spikes values
    spike = 2 * pi;
    YY = mod(Y, 2*pi);
    A1 = find(abs(YY(:, 1) - spike) < error);
    A2 = find(abs(YY(:, 2) - spike) < error);
    A1 = Find_Near_Points(A1);
    A2 = Find_Near_Points(A2);
    
    % Death condition
    if ((length(A1) < 1) && (length(A2) < 1))
        err = 3; % Both dead
    elseif (length(A2) < 1)
        err = 2; % Only second dead
    elseif (length(A1) < 1)
        err = 1; % Only first dead
    else
        err = 0; % Spiking process
    end

end

function AA1 = Find_Near_Points(A1)
    AA1 = A1;
    tolerance = 2;
    diffs = diff(A1);
    indices_to_remove = find(diffs < tolerance);
    AA1(indices_to_remove + 1) = [];
end

function AA1 = id(AA1, T)
    AA1 = num2cell(AA1);
    AA1 = T([AA1{:}], :);
end