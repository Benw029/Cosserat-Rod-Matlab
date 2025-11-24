function Tt = SMCCosserat(P, Pd, Vt, Vtd, C, k, epsilon, R, q, Ad)

    X1 = P;
    X1dot = R*q;

    
    e = Pd - P;
    edot = Vtd - Vt; %Xdotd - X1dot
    S = edot + C*e;

    alpha = 

    ac = (1/rho*A)*(n + fe);
    bc = (-1/rho*A)*(alpha);

    Tt = (C * (Vtd - Vt) + Ad - ac + epsilon*sign(S) + k*S)/bc;
end