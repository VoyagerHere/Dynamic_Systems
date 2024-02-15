% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%  Freq synchronization modes for phi-alpha system for 2-neurons
%
%  alpha - offset parameter of system
%  d - sync parameter of system
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

format longG
n = [3; 3; 3];

delta_max = 0.01;

g1 = 1.01;
g = [g1; g1 + delta_max; g1 + 2 * delta_max];

alpha = pi/8;
color = 'green';
global alpha_text
alpha_text = '2π/3';

alpha = 0;
color = 'green';
global alpha_text
alpha_text = '0';

DiagOut = diag(alpha, g, n);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function DiagOut = diag(alpha, g, n)
    d_list = 0:0.00001:0.02;
    DiagOut = zeros(length(d_list), 5);

   parfor m = 1:length(d_list)
        d = [d_list(m); d_list(m)];
        % DEBUG
        % disp(d);
        y0 = zeros(size(n, 1), 1);
        a = 4000;
        b = 10000;
        opts = odeset('RelTol', 2e-13, 'AbsTol', 1e-100, 'Refine', 10);
        [T, Y] = ode45(@(t, y)eqn(t, y, g, n, d, alpha), [a, b], y0, opts);
% draw(T,Y,g,d);
        [freq_ratio, err] = FREQ(T, Y, n);
        DiagOut(m, :) = [freq_ratio(1), freq_ratio(2), d(1), d(2), err];
   end
end

function [ratio, err] = FREQ(T, Y, n)
    TT = ms(Y,T,n);
    TT1 = nonzeros(TT(:,1));
    TT2 = nonzeros(TT(:,2));
    TT3 = nonzeros(TT(:,3));

    err = 0;

    if ((length(TT1) < 1) && (length(TT2) < 1) && (length(TT3) < 1))
        err = 4;
    elseif (length(TT3) < 1)
        err = 3;
    elseif (length(TT2) < 1)
        err = 2;
    elseif (length(TT1) < 1)
        err = 1;
    end
 
    if (err > 0)
        ratio = [0; 0];
        return;
    end

%     TT1 = del_middle_burst_int(n(1), TT1);
%     TT2 = del_middle_burst_int(n(2), TT2);
%     TT3 = del_middle_burst_int(n(2), TT2);
%     TT1 = mean(TT1);
%     TT2 = mean(TT2);
%     TT3 = mean(TT3);

    TT1 = sum(TT1);
    TT2 = sum(TT2);
    TT3 = sum(TT3);
 
    w1 = vpa(2 * pi / TT1);
    w2 = vpa(2 * pi / TT2);
    w3 = vpa(2 * pi / TT3);
    ratio =  [w1/w2; w2/w3];
 end

 function OUT = ms(Y,T,n)
  OUT = zeros(max(n),length(n));
  for ind = 1:length(n)
      UNSTBL = delete_unstb(T,Y, n, 10, ind);
      for k = 1:n(ind)
          [A,err] = find_spikes(Y, n, ind);
          if ((err) > 0)
              return;
          end
          start = find(abs((A - UNSTBL)) < 10);
          if (isnan(UNSTBL))
              start = 1;
          end
          Diff_points = A(start:end, 1);
          p1 = Diff_points(k);
          p2 = Diff_points(k+1);
          OUT(k, ind) = id(p2,T) - id(p1,T);
      end
  end
end

 function UNSTBL = delete_unstb(T, Y, n, error, ind)
  UNSTBL = NaN;
  [A, err] = find_spikes(Y, n, ind);
  if ((error == 0) && (err > 0))
      return;
  end
  num = n(ind);
  B = zeros(floor(size(A,1)/num),1);
  for m = 1:num:length(A)
      B(ceil(m/num),1) = A(m,1);
  end
  if ((length(B)) < (error + 2))
      return;
  end
  UNSTBL = B(error, 1);
end

function T_DIFF = del_middle_burst_int(n, T)
    if (n > 1)
        T = sort(T, 'descend');
        T_DIFF = T(2:size(T));
    else
        T_DIFF = T;
    end
end

% NAMED ms previously
function OUT = FIND_DIFF(Y, T, n, ind)
    tilda = tld(n);
    OUT = [];
    YY = mod(Y, 2 * pi);

    for k = 1:n
        Diff_points = Find_Near_Points(find(abs(YY(:, ind) - tilda) < 0.05));

        if (length(Diff_points) < n + 1)
            return;
        end

        A1 = Diff_points(k);
        A2 = Diff_points(k + 1);
        AA1 = id(A1, T);
        AA2 = id(A2, T);
        OUT(end + 1, :) = AA2 - AA1;
    end

end

function [A, err] = find_spikes(Y, n, ind)
  tilda = tld(n(ind));
  error = 0.05;
  YY = mod(Y, 2 * pi);
  A = find(abs(YY(:, ind) - tilda) < error);
  A = Find_Near_Points(A);

  if (find(abs(YY(:, ind) - 2*pi) < error))
    err = 0;
  else
    err = 1;
  end
end

function tilda = tld(n)
  tilda = n*pi/2 - floor(n/4)*2*pi;
  if abs(tilda - pi/2) < 0.02
      tilda = 3*pi/2;
  else
      tilda = abs(tilda - pi);
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

function dy_dt = eqn(t, y, g, n, d, alpha)
    f = g - sin(y ./ n);
    exch = d(1) * sin(y([2; 3; 1]) - y([1; 2; 3]) - alpha) + d(2) * sin(y([3; 1; 2]) - y([1; 2; 3]) - alpha);
    dy_dt = f + exch;
end

function draw(T, Y, g, d)
    global alpha_text;
    delta_max = g(2) - g(1);
    Y = mod(Y, 2 * pi);
    plot(T, Y(:, 1), '-b', T, Y(:, 2), '-g', T, Y(:, 3),'-r');
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
         sprintf('\\gamma_3 = %.2f', g(3)), ...
        'Interpreter', 'tex');
end
