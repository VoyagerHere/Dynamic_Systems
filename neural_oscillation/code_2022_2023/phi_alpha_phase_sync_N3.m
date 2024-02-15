% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%  Phase synchronization modes for phi-alpha system
%  
%  N = 3 - number of elements in system
%  alpha - offset parameter of system
%  d1, d2 - sync parameter of system
%
%  If n = 1 then burst is equivalent to spike
%
%  System state code declaration:
%       0** - error code  (death some of elements in system)
%       1** - global sync state (no sync, burst or spike sync)
%       2** - local sync state (between some of elements in system)
% 
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

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
% Parametrs to set
error = 30;
n1 = 3;
n2 = 3;
n3 = 3;
n = [n1; n2; n3];
g_num = 500;


global alpha_text
alpha_text = 'π/8';
alpha_text_par = 'pi_8';
alpha = pi / 8;
d_max = 0.05;  d_accuracy = 0.0001;


% alpha_text = '2π/3';
% alpha_text_par = '2pi_3';
% alpha = 2*pi / 3;
% d_max = 0.04;  d_accuracy = 0.0001;

% alpha_text = '0';
% alpha_text_par = '0';
% alpha = 0;
% d_max = 0.15;  d_accuracy = 0.0001;

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Initial parametrs
delta_phi_error = 0.05;
g1 = 1.01;
global delta_max
delta_max = 0.025;

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Main program
g_list = linspace(1.01, g1 + delta_max, g_num);
GAMMA_TABLE = [];
GAMMA_TABLE = gamma_iterate(g1, n, g_list, GAMMA_TABLE, ...
    alpha, d_max, d_accuracy, error);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Save file
filename = 'alpha_phase_sync';
folder = 'plots_data_131223';
alpha_text_par = sprintf('%s_%d_%d_%d', alpha_text_par, n1, n2, n3);
filename = sprintf('%s_%s', filename, alpha_text_par);
time = datestr(clock,'YYYY_mm_dd_HH_MM_SS');
filename = sprintf('%s/%s_%s.mat', folder, filename, time);
save(filename, '-v7.3', '-nocompression')


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Finctions
function GAMMA_TABLE = gamma_iterate(g1, n, g_list, ...
    GAMMA_TABLE, alpha, d_max, d_accuracy, error)
    d_list = 0:d_accuracy:d_max;
    for k = 1:length(g_list)
        g2 = g_list(k);
        delta = g2 - g1;
        g3 = g2 + delta;
        D_TABLE = d_iterate([g1; g2; g3], d_list, n, ...
            alpha, error);
        GAMMA_TABLE(end + 1:end + size(D_TABLE, 1), 1:7) = D_TABLE;
        disp(k);
    end

end

function DELTA_TABLE = d_iterate(g, dd, n, alpha, error)
    DELTA_TABLE = zeros(length(dd), 7);
    parfor m = 1:length(dd)
        d1 = dd(m);
        d2 = d1;
        y0 = [0.0; 0.0; 0.0];
        a = 8000;
        b = 20000;
        [T, Y] = ode45(@(t, y)eqn(t, y, g, n, d1, d2, alpha), [a, b], y0, opts);
% draw(T, Y, g, d1);
        [ratio1, ratio2, err, sync] = phase_sync(T, Y, n, error);
        delta = g(2) - g(1);
        DELTA_TABLE(m, :) = [d1, d2, ratio1, ratio2, delta, sync, err];
    end
end

function [ratio1, ratio2, err, sync] = phase_sync(T, Y, n, error)
    delta_phi_error = 0.005;
    [A1, err1] = find_spikes(Y, n, 1);
    [A2, err2] = find_spikes(Y, n, 2);
    [A3, err3] = find_spikes(Y, n, 3);
    err = detect_death([err1, err2, err3]);

    if (err == 111) % death
        DIFF_SP = NaN;
        DIFF_BS = NaN;
        ratio1 = -err;
        ratio2 = -err;
        sync = 100;
        return;
    end

    UNSTBL1 = delete_unstb(A1, err1, n, 1, error);
    UNSTBL2 = delete_unstb(A2, err2, n, 2, error);
    UNSTBL3 = delete_unstb(A3, err3, n, 3, error);

    A1 = A1(UNSTBL1:end);
    A2 = A2(UNSTBL2:end);
    A3 = A3(UNSTBL3:end);

    B1 = find_burst(A1, n(1));
    B2 = find_burst(A2, n(2));
    B3 = find_burst(A3, n(2));

    ratio1 = find_ratio(B1, B2);
    ratio2 = find_ratio(B2, B3);


    DIFF_BS_12 = find_diff(B1, B2, T);
    DIFF_BS_23 = find_diff(B2, B3, T);
    DIFF_BS = [DIFF_BS_12; DIFF_BS_23];

    DIFF_SP_12 = find_diff(A1, A2, T);
    DIFF_SP_23 = find_diff(A2, A3, T);
    DIFF_SP = [DIFF_SP_12; DIFF_SP_23];

    % sync code agreement
    % 2AB: 
    % A - burst(0) of spike(1) sync
    % B - configuration of synced elements
    sync = 100;
    if (any(el_delta(abs(DIFF_BS), delta_phi_error)) && any(err  == 0)) 
     sync = 101;
    elseif (el_delta(abs(DIFF_BS_12), delta_phi_error))
     sync = 202;
    elseif (el_delta(abs(DIFF_BS_23), delta_phi_error))
     sync = 201;
    end
    
    if (any(el_delta(abs(DIFF_SP), delta_phi_error)) && any(err  == 0))
     sync = 111;
    elseif (el_delta(abs(DIFF_SP_12), delta_phi_error))
     sync = 212;
    elseif (el_delta(abs(DIFF_SP_23), delta_phi_error))
     sync = 211;
    end
end

function B = find_burst(A, n)
    B = zeros(floor(size(A, 1) / n), 1);

    for m = n:n:length(A)
        B(round(m / n), 1) = A(m, 1);
    end

end

function out = el_delta(DIFF, error)
    mn = mean(abs(DIFF));
    out = abs(abs(DIFF) - mn) < error;
end

function ratio = find_ratio(A, B)
    ratio = zeros(length(A(:, 1)) - 2, 1);

    for i = 2:length(A) - 1
        ratio(i - 1, 1) = sum(B < A(i + 1, 1)) - sum(B < A(i, 1));
    end
    ratio = mode(ratio);
    if (ratio == 0)
        ratio = 1;
    end
end

function DIFF = find_diff(Arr1, Arr2, T)
    NEAR = []; % find near spikes

    for i = 1:length(Arr1)
        [~, index] = max(Arr2(Arr2 <= Arr1(i)));

        if ~isempty(index)
            NEAR(end + 1) = Arr2(index);
        end

    end

    Arr1 = id(Arr1, T);
    NEAR = id(NEAR, T);
    len = length(NEAR(:, 1));
    DIFF = zeros(len, 1);
    Arr1 = Arr1((end - len + 1):end, 1);
    DIFF(:, 1) = Arr1 - NEAR;
end

function [A, err] = find_spikes(Y, n, ind)
    tilda = tld(n(ind));
    error = 0.05;
    YY = mod(Y, 2 * pi);
    A = find(abs(YY(:, ind) - tilda) < error);
    A = Find_Near_Points(A);

    if (find(abs(YY(:, ind) - 2 * pi) < error))
        err = 0;
    else
        err = 1;
    end
end

function err = detect_death(errors)

    if ((errors(1)) && (errors(2)) && (errors(3)))
        err = 0111;
    elseif (errors(3))
        err = 0111;
    elseif (errors(2))
        err = 0110;
    elseif (errors(1))
        err = 0100;
    else
        err = 0;
    end
end

function UNSTBL = delete_unstb(A, err, n, ind, error)
    UNSTBL = 1;

    if ((error == 0) && (err > 0))
        return;
    end

    num = n(ind);
    B = zeros(floor(size(A, 1) / num), 1);

    for m = 1:num:length(A)
        B(ceil(m / num), 1) = A(m, 1);
    end

    if ((length(B)) < (error + 2))
        return;
    end

    unstbl_el = B(error, 1);
    UNSTBL = find(A == unstbl_el);

    if (isnan(UNSTBL))
        UNSTBL = 1;
    end

end

function draw(T, Y, g, d)
    global delta_max;
    global alpha_text;
    Y = mod(Y, 2 * pi);
    plot(T, Y(:, 1), '-b', T, Y(:, 2), '-g', T, Y(:, 3), '-r');
    ylim([0 2 * pi]); yticks([0 pi 2 * pi]); yticklabels(["0" "\pi" "2\pi"]);
    xlim([T(1, 1) T(end, 1)])
    xticks(linspace(T(1, 1), T(end, 1), 10)); % set to show only int
    ax = gca;
    ax.XTick = unique(round(ax.XTick));
    yticks([0 pi 2 * pi]);
    xlabel("t");
    % xlim([556 T(end,1)])
    ylabel('\phi', 'Interpreter', 'tex');
    title(sprintf('d = %g, Δ = %g, α = %s', d, delta_max, alpha_text));
    legend(sprintf('\\gamma_1 = %.2f', g(1)), ...
        sprintf('\\gamma_2 = %.2f', g(2)), ...
        'Interpreter', 'tex');
end

function dy_dt = eqn(t, y, g, n, d1, d2, alpha)
  f = g - sin(y ./ n);
  exch = [d1 * sin(y(2) - y(1) - alpha) + d2 * sin(y(3) - y(1) - alpha); d1 * sin(y(1) - y(2) - alpha) + d2 * sin(y(3) - y(2) - alpha); d1 * sin(y(1) - y(3) - alpha) + d2 * sin(y(2) - y(3) - alpha)];
  dy_dt = f + exch;
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

function tilda = tld(n)
    tilda = n * pi / 2 - floor(n / 4) * 2 * pi;

    if abs(tilda - pi / 2) < 0.02
        tilda = 3 * pi / 2;
    else
        tilda = abs(tilda - pi);
    end

end
