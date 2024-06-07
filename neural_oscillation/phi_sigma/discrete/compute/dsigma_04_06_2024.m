core = feature('numcores');
pool = parpool('local',core);
disp(['Pool has been started with Num Workers ' num2str(pool.NumWorkers)]);

retries = 0;
retry_limit = 3;
while (pool.NumWorkers < core)
    retries = retries + 1;
    disp('Restarting parallel pool');
    delete(pool);
    pool = parpool('local',core);
    disp(['Pool has been started with Num Workers ' num2str(pool.NumWorkers)]);
    if(retries >= retry_limit)
        break;
    end
end


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% MAIN PROGRAM
% 
n1 = 1;
n2 = 1;
n = [n1; n2];

delta_phi_error = 0.05;
DELTA = 0.01;
g1 = 1.001;
g2 = g1 + DELTA;

disp("Start calculation");

% Table with gamma
sigma_list = linspace(0, 1, 400);
SIGMA_TABLE = [];
SIGMA_TABLE = sigma_iterate(g1, g2, n, sigma_list, SIGMA_TABLE);
close all;

% Save file
time = datestr(clock,'YYYY_mm_dd_HH_MM_SS');
filename = 'dsigma';
filename = sprintf('%s_%s.mat', filename, time);
save(filename, '-v7.3', '-nocompression')


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% FUNCTIONS
% 

% Iterate over the sigma par value
function SIGMA_TABLE = sigma_iterate(g1, g2, n, sigma_list, SIGMA_TABLE)
    delta_phi_error = 0.05; % phase sync tolerance
    d_list = 0:0.0005:0.1;
    num_of_iterations = length(sigma_list); 
    close all;
    for k = 1:length(sigma_list)
        sigma = sigma_list(k);
        D_TABLE = d_iterate([g1; g2], d_list, n, sigma, delta_phi_error);
        SIGMA_TABLE(end + 1:end + size(D_TABLE, 1), 1:3) = D_TABLE;
        fprintf('Iteration = %d, All = %d \n',k, num_of_iterations);
        drawnow('update');
        disp(k);
    end
end

% Iterate over d values
function DELTA_TABLE = d_iterate(g, d_list, n, sigma, delta_phi_error)
    DELTA_TABLE = zeros(length(d_list), 3);
    
   parfor m = 1:length(d_list)
        d = d_list(m);
        y0 = [pi/2; pi/2];
        a = 8000;
        b = 10000;
        step = 0.001;
        % Integrate the system with one-step Euler method
        %[T, Y] = Euler_int_syst(d, sigma, n, y0, a, b, step);
        F1_0 = 0;
        F2_0 = 0;
        F = [F1_0; F2_0];
        [T, Y] = odeK4(@eqn ,[a,b-step], y0, step,[42, 1] ,{d, n, sigma, F, g}, 0, 0);
        T = transpose(T);
        A1 = Y(:,1); % Phi_1 values
        A2 = Y(:,2); % Phi_2 values
        [A1, A2, err] = find_spikes(Y);

        %if (err == 0)
            %[T, Y] = delete_unstb(T, Y, A1, A2, 6);
        %end
        [DIFF, diff, ratio] = ms_phase_sync(T, Y);
% DEBUG   
% draw(T,Y, diff, g, d, sigma);
        if (el_delta(abs(DIFF), delta_phi_error)) % phase sync
            DELTA_TABLE(m, :) = [sigma, d, ratio];
        elseif (err > 0) % death
            DELTA_TABLE(m, :) = [sigma, d, ratio];
        else %  quasi-periodic
            DELTA_TABLE(m, :) = [sigma, d, -9];
        end
    end
end

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
  xout(ii,:) = transpose(x);
  %xdot(:,ii) = xdotRK4;
  %%Integrate
  
  if extraparameters
    if find(y == 42)
      [xdot1,outs] = feval(functionHandle, time, x, parameters{1}, parameters{2}, parameters{3}, parameters{4}, parameters{5});
      other(:,ii) = outs;
      newF = outs;
      [xdot2,outs] = feval(functionHandle, time + (.5*timestep), x + (xdot1*.5*timestep), parameters{1}, parameters{2}, parameters{3}, parameters{4}, parameters{5});
      [xdot3,outs] = feval(functionHandle, time + (.5*timestep), x + (xdot2*.5*timestep), parameters{1}, parameters{2}, parameters{3}, parameters{4}, parameters{5});
      [xdot4,outs] = feval(functionHandle, time + timestep, x + (xdot3*timestep), parameters{1}, parameters{2}, parameters{3}, parameters{4}, parameters{5});
      
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

% Integrate the system with one-step Euler method
function [T, Y] = Euler_int_syst(d, sigma, n, y0, a, b, step)
  h = 1/step;
  T = a:h:(b-h);
  T = transpose(T);
  Y(length(T), 2) = 0;
  Y(1,1:2) = y0;
  F1_0 = 0;
  F2_0 = 0;
  F = [F1_0; F2_0];

  % Solve ODE for each time step and update F
  for i = 1:((length(T)-1))
      [dy_dt, F] = eqn(T(i), Y(i,:), d, n, sigma, F);
      Y(i+1,1) = Y(i,1) + dy_dt(1)*(h);
      Y(i+1,2) = Y(i,2) + dy_dt(2)*(h);
  end
end

% Find spikes on phi values
function [A1, A2, err] = find_spikes(Y)
    error = 0.05; % Tolerance for spikes values
    spike = 2 * pi;
    YY = mod(Y, 2 * pi);
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

% Increse integration right interval for
% more correct result on death
function b = set_dynamic_death(g, d, def_b)
    if ((g(2) < 0.006) && (d > 0.025))
        b = 10000;
    elseif ((g(2) > 0.006) && (d > 0.036))
        b = 10000;
    else
        b = def_b;
    end
end

% Delete spikes on transit process
function [T, Y] = delete_unstb(T, Y, A1, A2, error)
    % error - count of bad spike
    % default  error = 5;
    %   draw(T,Y,0);
    error = min(error, max((length(A1) - 4), 0));
    if (error == 0)
        return;
    end
    first_spike = A1(error, 1);
    Y = Y(first_spike:end, :);
    T = T(first_spike:end, :);
end

% Check if value correct with error
function out = el_delta(DIFF, error)
    mn = mean(abs(DIFF));
    out = abs(abs(DIFF) - mn) < error;
end

% Find spikes ratio
% DIFF values - values between nearest spikes to left
% RARIO values - number of spikes of second el between two spikes of
%   first element
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

% Phase realisation function (for DEBUG)
function draw(T, Y, diff, g, d, sigma)
    Y = mod(Y, 2 * pi);
    plot(T, Y(:, 1), '-r', T, Y(:, 2), '-b');
    ylim([0 2 * pi]); yticks([0 pi 2 * pi]); yticklabels(["0" "\pi" "2\pi"]);
    xlim([T(1, 1) T(end, 1)])
    xticks(linspace(T(1, 1), T(end, 1), 10)); % set to show only int
    ax = gca;
    ax.XTick = unique(round(ax.XTick));
    yticks([0 pi 2 * pi]);
    xlabel("t");
    % xlim([556 T(end,1)])
    ylabel('\phi', 'Interpreter', 'tex');
    title(sprintf('d = %g, σ = %g, |φ₁ - φ₂| ≈ %.2f', d, sigma, diff));
    legend(sprintf('\\gamma_1 = %g', g(1)), ...
        sprintf('\\gamma_2 = %g', g(2)), ...
        'Interpreter', 'tex');
end


function [dy_dt, F] = eqn(~, y, d, no, sigma, F, g)
  yy = mod(y, 2*pi);
  
  f = g - sin(y ./ no);
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