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



%  addpath('additional\progress_bar\')
% addpath('additional\')


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Parametrs to set
error = 20;
n1 = 2;
n2 = 2;
n = [n1; n2];
g_num = 500;


% global alpha_text
% alpha_text = 'π/8';
% alpha_text_par = 'pi_8';
% alpha = pi / 8;
% d_max = 0.07;  d_accuracy = 0.0001;


alpha_text = '2π/3';
alpha_text_par = '2pi_3';
alpha = 2*pi / 3;
d_max = 0.09;  d_accuracy = 0.00001;


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
alpha_text_par = sprintf('%s_%d_%d_%d', alpha_text_par, n1, n2);
filename = sprintf('%s_%s', filename, alpha_text_par);
time = datestr(clock,'YYYY_mm_dd_HH_MM_SS');
filename = sprintf('%s_%s.mat', filename, time);
save(filename, '-v7.3', '-nocompression')



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Finctions
function GAMMA_TABLE = gamma_iterate(g1, n, g_list, ...
    GAMMA_TABLE, alpha, d_max, d_accuracy, error)
    delta_phi_error = 0.005;
    d_list = 0:d_accuracy:d_max;
    num_of_iterations = length(g_list);
    for k = 1:length(g_list)
        g2 = g_list(k);
        D_TABLE = d_iterate([g1; g2], d_list, n, ...
            alpha, delta_phi_error, error);
        GAMMA_TABLE(end + 1:end + size(D_TABLE, 1), 1:4) = D_TABLE;
        disp(num_of_iterations);
    end

end

function DELTA_TABLE = d_iterate(g, dd, n, alpha, delta_phi_error, error)
    DELTA_TABLE = zeros(length(dd), 4);
    parfor m = 1:length(dd)
        d = dd(m);
        y0 = [0; 0];
        a = 2000;
        b = 6000;
        opts = odeset('RelTol',1e-13,'AbsTol',1e-14,'Refine',10);
        [T, Y] = ode45(@(t, y)eqn(g, n, t, y, d, alpha), [a, b], y0, opts);
% draw(T, Y, g, d);
        [DIFF_SP, DIFF_BS, ratio, err] = phase_sync(T, Y, n, error);
        delta = g(2) - g(1);

        sync_mode = 100; % global sync state % no sync

        if (el_delta(abs(DIFF_BS), delta_phi_error)) % burst sync
            sync_mode = 101;
        end

        if (el_delta(abs(DIFF_SP), delta_phi_error)) % spike sync
            sync_mode = 111;
        end

        if (sync_mode > 100) % sync
            DELTA_TABLE(m, :) = [d, ratio, delta, sync_mode];
        elseif (err > 0) % death
            DELTA_TABLE(m, :) = [d, ratio, delta, sync_mode];
        else %  quasi-periodic
            DELTA_TABLE(m, :) = [d, 0, delta, sync_mode];
        end

    end

end

function [DIFF_SP, DIFF_BS, ratio, err] = phase_sync(T, Y, n, error)
    [A1, err1] = find_spikes(Y, n, 1);
    [A2, err2] = find_spikes(Y, n, 2);
    err = detect_death([err1, err2]);

    if (err > 000) % death
        DIFF_SP = NaN;
        DIFF_BS = NaN;
        ratio = -err;
        return;
    end

    UNSTBL1 = delete_unstb(A1, err1, n, 1, error);
    UNSTBL2 = delete_unstb(A2, err2, n, 2, error);
    A1 = A1(UNSTBL1:end);
    A2 = A2(UNSTBL2:end);
    B1 = find_burst(A1, n(1));
    B2 = find_burst(A2, n(2));

    ratio = find_ratio(B1, B2);
    DIFF_BS = find_diff(B1, B2, T);
    DIFF_SP = find_diff(A1, A2, T);
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

    if (errors(1)) && (errors(2))
        err = 11;
    elseif (errors(2))
        err = 10;
    elseif (errors(1))
        err = 01;
    else
        err = 0;
    end
end

function UNSTBL = delete_unstb(A, err, n, ind, error)
    UNSTBL = 1;

    if ((error == 0) || (err > 0))
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
    plot(T, Y(:, 1), '-b', T, Y(:, 2), '-g');
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

function dy_dt = eqn(g, n, ~, y, d, alpha)
    f = g - sin(y ./ n);
    exch = [d * sin(y(2) - y(1) - alpha); d * (sin(y(1) - y(2) - alpha))];
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
