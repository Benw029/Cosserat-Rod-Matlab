% This program implements the Cosserat Rod Theory for a Tendon Driven CR
% with simple Control inputs using Sliding Model Control and Lyapunov
% Stability
% Benjamin Chung
% 11/25/2025

function SingleTendonCosseratRod
clear all;
clc;

global i p R j n m v u q w ns vs us vt ut qt wt vst ust vh uh vsh ush qh wh nLL mLL x y z X Y Z  %Make vars available in whole program
global Tension_Hist X1_Hist Xdot1_Hist e_Hist edot_Hist f_Hist norm_e_Hist Control_c
%Hat operator
hat = @(y)[0, -y(3), y(2);
    y(3), 0, -y(1);
    -y(2), y(1), 0];

%Declare variables
L = 0.5;                       %Length in m
N = 20;                        %Spatial resolution
E = 190e9;                     %Young's modulus
r = 1/1000;                     %Cross-section radius
rt = {[0.01;0;0]               %Location of Tendons - We'll see if 4 works
    [0;0.01;0]
    [-0.01;0;0]
    [0;-0.01;0]};
rho = 17189;                    %Density
g = [0;0;-9.81];               %Gravity vector
Bse = zeros(3);                %Shear/Bending Coefficient
Bbt = 0.008*eye(3);              %Bending/Torsion Coefficient
C = 0.01*eye(3);                %Viscous Damping Coefficient
Tt = {0 0 0 0};                %Assume 0 Tension during initial Static solve
dt = 0.005;                    %Size of Time step
alpha = -0.2;                  %BDF-alpha parameter
STEPS = 100;                    %Number of timesteps to completion
STEPCOUNT = 1;

%Initial Pre-bending variables
vstar = @(s)[0;0;1];
ustar = @(s)[0;0;0];
vsstar = @(s)[0;0;1];
usstar = @(s)[0;0;0];

%Boundary Conditions
%Clamped base to free end
for i = 1 : STEPS
    p{i,1} = [0;0;0];
    R{i,1} = eye(3);
    q{i,1} = [0;0;0];
    w{i,1} = [0;0;0];
end

nL = [0;0;0];                    %Start with no forces
mL = [0;0;0];

%Dependent Parameter Calculations
A = pi*r^2;                                 %Cross-sectional area
J = diag([pi*r^4/4  pi*r^4/4  pi*r^4/2]);   %Inertia Matrix
G = E / (2 * (1 + 0.3));                    %Shear modulus
Kse = diag([G*A, G*A, E*A]);                %SE Stiffness Matrix
Kbt = diag([E*J(1,1), E*J(2,2), G*J(3,3)]); %BT Stiffness Matrix
ds = L/(N-1);                               %Finite change in S over total steps

%BDF-alpha coefficients
C0 = (1.5 + alpha) / ( dt*(1+alpha));
C1 = -2/dt;
C2 = (0.5 + alpha) / ( dt*(1+alpha));
d1 = alpha / (1+alpha);

% -- START SIMULATION --
% Steps to consider:
% 1) Solve with initial static conditions
% 2) Apply Tension Forces
% 3) Semi-discretize the PDE into an ODE using BDF-Alpha
% 4) Solve PDE with shooting

%Initial Static Solve in order to get IVP problem going
i = 1;
disp("Solving Initial Static Model")
fsolve(@initialSolve, zeros(6,1)); %Solve static BVP w/ shooting method
visualize();

%Setup initial history terms
for j = 1 : N-1
    vh{i+1,j} = (C1+C2)*v{i,j};
    uh{i+1,j} = (C1+C2)*u{i,j};
    vsh{i+1,j} = (C1+C2)*vs{i,j};
    ush{i+1,j} = (C1+C2)*us{i,j};
    qh{i+1,j} = [0;0;0];
    wh{i+1,j} = [0;0;0];
    q{i,j} = [0;0;0];
    w{i,j} = [0;0;0];
end

%Set Control Variables
CONTROL = true;                     %SET TRUE OR FALSE IF YOU WANT TO USE THE SMC
P_d = [0.027; 0.0000; 0.4727];      %Desired end point
V_d = [0; 0; 0];                    %Desired tip velocity
Acc_d = [0; 0; 0];                  %Desired tip acceleration
Control_c = diag([15000, 1, 500]);  %Control gain 1
k = diag([60, 1, 5]);               %Control gain 2 (for S convergence)
epsilon = 0.005;                    %Reaching epsilon
clamp = 8;

%DEBUG
warning('off', 'all');
options = optimoptions('fsolve', 'Display', 'none');


% -- START SIMULATION LOOP --
for i = 2 : STEPS
    STEPCOUNT = STEPCOUNT + 1;
    %Set the Tendon Forces at each time step
    if i < 5
        Tt{1} = 0; %Apply 0 for initial history terms
        Tt{2} = 0;
        Tt{3} = 0;
        Tt{4} = 0;
    else
        if CONTROL
            Tt{1} = SMCCosserat(P_d, V_d, Control_c, k, epsilon, Acc_d);
            Tt{2} = 0;
            Tt{3} = 0;
            Tt{4} = 0;
        else
            Tt{1} = 3; %MANUAL INPUT TENDON FORCES
            Tt{2} = 0;
            Tt{3} = 0;
            Tt{4} = 0;
        end
    end

    %Now to actually solve, given the tendon forces
    disp(["Solving Dynamic Model, Step = ", i])
    fsolve(@dynamicSolve, [n{i-1,1}; m{i-1,1}], options); %Solve semi-discretized PDE w/ shooting
    updateBDFalpha();
    visualize();
end

disp('Final Pos'); p{i,N-1}

% -- FUNCTIONS --
    function E = initialSolve(G)
        n{i,1} = G(1:3);
        m{i,1} = G(4:6);

        %Euler's method for IVP problem
        for j = 1 : N-1
            [ps, Rs, ns, ms, us{i,j}, vs{i,j} ,v{i,j}, u{i,j}] = formStaticCosseratODE(p{i,j},R{i,j},n{i,j},m{i,j});
            p{i,j+1} = p{i,j} + ds*ps;
            R{i,j+1} = R{i,j} + ds*Rs;
            n{i,j+1} = n{i,j} + ds*ns;
            m{i,j+1} = m{i,j} + ds*ms;
        end

        E = [ n{i,N} - nL;  m{i,N} - mL ];
    end

    function E = dynamicSolve(G)
        n{i,1} = G(1:3);
        m{i,1} = G(4:6);

        %Euler's method
        for j = 1 : N-1
            [ps, Rs, ns, ms, qs, ws, vs{i,j}, us{i,j},...
                v{i,j}, u{i,j}, vt{i,j}, ut{i,j},...
                qt{i,j}, wt{i,j},vst{i,j}, ust{i,j}] = formCosseratODE(p{i,j},R{i,j},n{i,j},m{i,j},q{i,j},w{i,j});

            p{i,j+1} = p{i,j} + ds*ps;
            R{i,j+1} = R{i,j} + ds*Rs;
            n{i,j+1} = n{i,j} + ds*ns;
            m{i,j+1} = m{i,j} + ds*ms;
            q{i,j+1} = q{i,j} + ds*qs;
            w{i,j+1} = w{i,j} + ds*ws;

        end

        E = [n{i,N} - nLL ;  m{i,N} - mLL];
    end


%This function will essentially build the ODE that will be put through
%fsolve dynamically over time
    function [ps,Rs,ns,ms,qs,ws,vs,us,v,u,vt,ut,qt,wt,vst,ust] = formCosseratODE(p, R, n, m, q, w, Tt)
        v = (Kse + C0*Bse)\(R'*n + Kse*vstar(ds*(j-1)) - Bse*vh{i,j}); %Use \ for inverse for better stability
        u = (Kbt + C0*Bbt)\(R'*m + Kbt*ustar(ds*(j-1)) - Bbt*uh{i,j});

        %Time derivatives - Replaced by BDFa terms
        vt = C0*v + vh{i,j};
        ut = C0*u + uh{i,j};
        qt = C0*q + qh{i,j};
        wt = C0*w + wh{i,j};

        %Calculate Forces
        %Using _ for subscript t, for tension forces
        [a_t, b_t, A_t, B_t, G_t, H_t] = getTendonForces(v, u);

        phi = [Kse + C0*Bse + A_t, G_t;
            G_t', Kbt + C0*Bbt + H_t];

        LamdaN = -a_t + rho*A*(hat(w)*q + qt)+ C*q.*abs(q)- R'*rho*A*g;
        LamdaM = -b_t + rho*(hat(w)*J*w + J*wt) -hat (v)*(Kse*(v-vstar(ds*(j-1)))+Bse*vt);
        GammaV = hat(u)*(Kse*(v-vstar(ds*(j-1)))+Bse*vt)-Kse*vsstar(ds*(j-1))+Bse*vsh{i,j};
        GammaU = hat(u)*(Kbt*(u-ustar(ds*(j-1)))+Bbt*ut)-Kbt*usstar(ds*(j-1))+Bbt*ush{i,j};

        %Matrix inverse formula as all matrices are orthogonal
        us = 1/det(phi)*(-G_t'*(-GammaV+LamdaN)+(Kse+C0*Bse+A_t)*(-GammaU+LamdaM));
        vs = 1/det(phi)*((Kbt+C0*Bbt+H_t)*(-GammaV+LamdaN)-G_t*(-GammaU+LamdaM));

        vst = C0*vs + vsh{i,j};
        ust = C0*us + ush{i,j};

        f = rho*A*g + R*(a_t + A_t*vs + G_t*us) - (R*C*q.*abs(q));
        f_Hist{i,j} = f; 
        l = R*(b_t + B_t*vt + H_t*ut);

        %Find boundary forces and moments from tendons for residual calculation
        [nLL, mLL] = getTendonForcesAtL(R, u, v);

        %Spatial derivatives
        ps = R * v;
        Rs = R * hat(u);
        qs = vt - hat(u)*q + hat(w)*v;
        ws = ut - hat(u)*w;
        ns = rho*A*R*(hat(w)*q + qt) - f;
        ms = R*rho*(hat(w)*J*w + J*wt) - hat(ps)*n - l;
    end

%STATIC
%This function is very explicit as I simly set the time terms to 0
%Could be optimized to OMIT time terms later on
    function [ps, Rs, ns, ms, us, vs, v, u] = formStaticCosseratODE(p, R, n, m)
        v = (Kse)\(R'*n + Kse*vstar(ds*(j-1))); %Use \ for inverse for better stability
        u = (Kbt)\(R'*m + Kbt*ustar(ds*(j-1)));

        %Time derivatives - all 0 at static

        %Calculate Forces
        %Using _ for subscript t, for tension forces
        [a_t, b_t, A_t, B_t, G_t, H_t] = getTendonForces(v, u);

        %Following Sharifi et al for vs and us
        %Their code seems to assume that all block matrices are diagonal
        %The following works in all cases, but is a bit less efficient
        %I didn't time to verify the matrix sizes to know if Gt is diagonal
        %WARNING: If singular, the whole solution explodes
        phi = [Kse + A_t, G_t;
               G_t', Kbt + H_t];

        LamdaN = hat(u) * (Kse * (v - vstar(ds*(j-1)))) - Kse*vsstar(ds*(j-1));
        LamdaM = hat(u) * (Kbt * (u - ustar(ds*(j-1)))) - Kbt*usstar(ds*(j-1));
        GammaN = R'*rho*A_t*g - a_t; %**R'R = 1, so can be omitted
        GammaM = hat(v) * Kse * (v - vstar(ds*(j-1))) - b_t; %**hat(ps) = Rv, also le = 0

        rhs = [-GammaN + LamdaN;
            -GammaM + LamdaM];
        temp = phi \ rhs;

        %Extract the result into vs and us
        vs = temp(1:3);
        us = temp(4:6);

        f = rho*A*g;
        f_Hist{i,j} = f;
        l = 0;

        %Spatial derivatives
        ps = R * v;
        Rs = R * hat(u);
        ns = -f;
        ms = -hat(ps)*n - l;
    end

%This function calculates the variables required to find the tendon forces
%Honestly, I could make this modular to have n tendons
%Probably a future thing
    function [a, b, At, Bt, Gt, Ht] = getTendonForces(v, u)

        %Tendon 1
        ptsb1 = hat(u)*rt{1}+v;
        At1 = -Tt{1} * (hat(ptsb1) * hat(ptsb1)) / (norm(ptsb1)^3);
        Gt1 = At1*hat(rt{1});
        Bt1 = hat(rt{1}) * At1;
        H1 = Bt1 * hat(rt{1});
        a1 = At1 * (hat(u) * ptsb1); %Assume the derivatives of r to be 0
        b1 = hat(rt{1}) * a1;

        %Tendon 2
        ptsb2 = hat(u)*rt{2}+v;
        At2 = -Tt{2} * (hat(ptsb2) * hat(ptsb2)) / (norm(ptsb2)^3);
        Gt2 = At2*hat(rt{2});
        Bt2 = hat(rt{2}) * At2;
        H2 = Bt2 * hat(rt{2});
        a2 = At2 * (hat(u) * ptsb2);
        b2 = hat(rt{2}) * a2;

        %Tendon 3
        ptsb3 = hat(u)*rt{3}+v;
        At3 = -Tt{3} * (hat(ptsb3) * hat(ptsb3)) / (norm(ptsb3)^3);
        Gt3 = At3*hat(rt{3});
        Bt3 = hat(rt{3}) * At3;
        H3 = Bt3 * hat(rt{3});
        a3 = At3 * (hat(u) * ptsb3);
        b3 = hat(rt{3}) * a3;

        %Tendon 4
        ptsb4 = hat(u)*rt{4}+v;
        At4 = -Tt{4} * (hat(ptsb4) * hat(ptsb4)) / (norm(ptsb4)^3);
        Gt4 = At4*hat(rt{4});
        Bt4 = hat(rt{4}) * At4;
        H4 = Bt4 * hat(rt{4});
        a4 = At4 * (hat(u) * ptsb4);
        b4 = hat(rt{4}) * a4;

        %Sum the values
        a = a1 + a2 + a3 +a4;
        b = b1 + b2 + b3 +b4;
        At = At1 + At2 + At3 + At4;
        Bt = Bt1 + Bt2 + Bt3 + Bt4;
        Gt = Gt1 + Gt2 + Gt3 + Gt4;
        Ht = H1 + H2 + H3 + H4;
    end

%This function find the tendon forces at the end point L
%Later then used in fsolve for residual calculation to solve the IVP
    function [nLL, mLL] = getTendonForcesAtL(R, u, v)
        pts1 = R*hat(u)*rt{1}+R*v;
        pts2 = R*hat(u)*rt{2}+R*v;
        pts3 = R*hat(u)*rt{3}+R*v;
        pts4 = R*hat(u)*rt{4}+R*v;

        nLL = -Tt{1}*pts1/norm(pts1) - Tt{2}*pts2/norm(pts2) - Tt{3}*pts3/norm(pts3); - Tt{4}*pts4/norm(pts4);
        mLL = -Tt{1}*hat(R*rt{1})*pts1/norm(pts1) - Tt{2}*hat(R*rt{2})*pts2/norm(pts2) - Tt{3}*hat(R*rt{3})*pts3/norm(pts3) - Tt{4}*hat(R*rt{4})*pts4/norm(pts4);
    end

    function updateBDFalpha()
        %h for history, these are history term updates
        for j = 1 : N-1
            vh{i+1,j} = C1*v{i,j} + C2*v{i-1,j} + d1*vt{i,j};
            uh{i+1,j} = C1*u{i,j} + C2*u{i-1,j} + d1*ut{i,j};
            vsh{i+1,j} = C1*vs{i,j} + C2*vs{i-1,j} + d1*vst{i,j};
            ush{i+1,j} = C1*us{i,j} + C2*us{i-1,j} + d1*ust{i,j};
            qh{i+1,j} = C1*q{i,j} + C2*q{i-1,j} + d1*qt{i,j};
            wh{i+1,j} = C1*w{i,j} + C2*w{i-1,j} + d1*wt{i,j};
        end
    end

    function Tt = SMCCosserat(X1d, Xdot1d, Control_c, k, epsilon, Xddotd)

        X1 = p{i-1,N-1};
        Xdot1 = R{i-1,N-1}*q{i-1,N-1};
        e = (X1d - X1);                 %Error between positions
        edot = Xdot1d - Xdot1;          %Derivative of error, which is velocity
        S = edot + Control_c*e;         %Sliding surface

        %Calculating the Alpha
        pts1 = R{i-1, N-1}*hat(u{i-1,N-1})*rt{1} + v{i-1,N-1};
        ptss1 = R{i-1, N-1}*hat(u{i-1,N-1})*hat(u{i-1,N-1})*rt{1} + hat(u{i-1,N-1})*v{i-1,N-1} + hat(us{i-1, N-1})*rt{1} + vs{i-1,N-1};

        a1 = ( (hat(pts1)*hat(pts1)*ptss1) / (norm(pts1)^3) ) + (pts1 / norm(pts1));
        alpha = a1;

        ac = (ns + rho*A*g - (R{i-1, N-1}*C*q{i-1, N-1}.*abs(q{i-1, N-1}))) / (rho*A);
        bc = (alpha) / (-rho*A);

        Tt = bc\(Control_c*edot + Xddotd - ac + epsilon*sign(S) + k*(edot + Control_c*e));

        %Apply clamping saturation
        if(Tt > clamp)
            Tt = clamp;
        end

        norm_e_Hist(i) = norm(e);
        Tension_Hist(i) = Tt;
        X1_Hist(i) = X1(1);
        Xdot1_Hist(i) = Xdot1(1);
        e_Hist(i) = e(1);
        edot_Hist(i) = edot(1);
    end


    function visualize()
        % p: cell array p{i,j} = 3x1 positions
        % N: number of points along rod
        % L: rod length for axis limits
        % i: current time step

        % Extract x,y,z positions
        x = zeros(1,N);
        y = zeros(1,N);
        z = zeros(1,N);
        for j = 1:N
            x(j) = p{i,j}(1);
            y(j) = p{i,j}(2);
            z(j) = p{i,j}(3);
        end

        % Persistent line handle to update instead of replotting
        persistent hLine
        if isempty(hLine) || ~isvalid(hLine)
            figure(1); clf;
            hLine = plot3(x, y, z, '-o', 'LineWidth', 2, 'MarkerSize', 4);
            axis([-0.1*L, 0.1*L, -0.1*L, 0.1*L, -0.05*L, 1.1*L]);
            xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');
            grid on; view(3);
        else
            % Update existing line
            set(hLine, 'XData', x, 'YData', y, 'ZData', z);
        end
        drawnow;
        pause(0.05);
    end

for time = 1:STEPCOUNT
    x(time) = p{time,N}(1);
    z(time) = p{time,N}(3);
    velocity{time} = R{time, N-1}*q{time,N-1};
end

velocity_x = cellfun(@(v) v(1), velocity);
velocity_y = cellfun(@(v) v(2), velocity);
velocity_z = cellfun(@(v) v(3), velocity);

accel_x = gradient(velocity_x, dt);
accel_y = gradient(velocity_y, dt);
accel_z = gradient(velocity_z, dt);



%All of the plotting fun!
if(CONTROL)
    x_axis = 1:size(x,2);
    figure (2)
    subplot(2,1,1)
    plot(x_axis,x, '-o');
    xlabel('Step');  ylabel('x (m)'); title('Tip Displacement - X Component');

    subplot(2,1,2)
    t_mean = mean(Tension_Hist);
    hold on
    plot_t = plot(x_axis, Tension_Hist, '-o');
    mean_t = yline(t_mean, 'r--', 'LineWidth', 2);
    xlabel('t (s)');  ylabel('F (N)'); title('Tendon Tension');
    legend([plot_t, mean_t], {'Tension (N)', ['Mean = ', num2str(t_mean,'%.2f')]}, 'Location','best')
    hold off

    figure(3)
    plot(X1_Hist, Xdot1_Hist);
    xlabel('X1'); ylabel('Xdot1');

    figure(4)
    plot(e_Hist, edot_Hist);
    xlabel('e'); ylabel('e_dot');

    figure(5)
    plot(x_axis, norm_e_Hist);
    xlabel('Step'); ylabel('e');

else
    x_axis = 1:size(x,2);
    figure (2)
    plot(x_axis, x, '-o');
    xlabel('Step');  ylabel('x (m)'); title('Tip Displacement - X Component');

    figure (3)
    plot(x_axis, z, '-o');
    xlabel('Step');  ylabel('z (m)'); title('Tip Displacement - Z Component');

    figure (4)
    hold on
    vx = plot(x_axis, velocity_x, 'r');
    vy = plot(x_axis, velocity_y, 'g');
    vz = plot(x_axis, velocity_z, 'b');
    xlabel('Step');  ylabel('Velocity (m/s)'); title('Tip Velocity');
    legend([vx vy vz], {'X Velocity','Y Velocity','Z Velocity'}, 'Location','best')
    hold off

    figure (5)
    hold on
    ax = plot(x_axis, accel_x, 'r');
    ay = plot(x_axis, accel_y, 'g');
    az = plot(x_axis, accel_z, 'b');
    xlabel('Step');  ylabel('Acceleration (m/s^2)'); title('Tip Acceleration');
    legend([ax ay az], {'X Acceleration','Y Acceleration','Z Acceleration'}, 'Location','best')
    hold off
end

end


