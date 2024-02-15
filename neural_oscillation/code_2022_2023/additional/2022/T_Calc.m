g1 = 1.02;
n = 5;
tilda = tld(n);

OUT = [];

for k = 1:n 
    phi1 = tilda+2*pi*(k-1);
    phi2 = phi1 + 2*pi;
    
    x1 =  (2 * n * (atan((g1*(tan(phi1/(2*n))- 1)/sqrt(g1^2 - 1)))))/sqrt(g1^2 - 1);
    x2 =  (2 * n * (atan((g1*(tan(phi2/(2*n))- 1)/sqrt(g1^2 - 1)))))/sqrt(g1^2 - 1);
    T = x2 - x1;
    x3 =  (2 * n * (atan((g1*(tan(phi2/(2*n))- 1)/sqrt(g1^2 - 1)))+pi))/sqrt(g1^2 - 1);
    TT = x3 - x1;
    if (T > 0)
      OUT(end+1, 1) = T;
    else
      OUT(end+1, 1) = TT;
    end
end

w = (2*pi) / sum(OUT);
% 
% 
% time1 = phi1 / pi;
% time2 = phi2 / pi;

function tilda = tld(n)
    tilda = n*pi/2 - floor(n/4)*2*pi;
    if abs(tilda - pi/2) < 0.02
        tilda = 3*pi/2;
    else
        tilda = abs(tilda - pi);
    end
end