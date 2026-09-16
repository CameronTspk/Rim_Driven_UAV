%% twin_vtol_params.m
%  Parameter set for the Twin Rim-Driven Vectored-Thrust UAV.
%  Run this before opening/simulating the Simulink model. All SI units.
%  Values marked [A]ssumed, [E]stimated, [M]easurable are defined in the
%  accompanying design document's Assumption Register -- replace [M] with
%  CAD- or bench-derived numbers as they become available.

p = struct();

% --- environment ---
p.g      = 9.81;              % [A] gravity                         [m/s^2]
p.rho    = 1.225;             % [A] air density (sea level)         [kg/m^3]

% --- mass / inertia ---  (replace with CAD)
p.m      = 1.50;             % [E] vehicle mass                     [kg]
p.I      = diag([0.07 0.02 0.075]);  % [E] inertia Ixx Iyy Izz      [kg m^2]
% NOTE Iyy << Ixx,Izz because both propulsors lie on the y (roll/yaw) axis.

% --- geometry (lateral baseline) ---
p.d      = 0.30;             % [M] half-span, propulsor +/-y offset [m]
p.r1     = [0;  p.d; 0];     %     right propulsor position (CG-rel)[m]
p.r2     = [0; -p.d; 0];     %     left  propulsor position (CG-rel)[m]

% --- propeller / duct ---
p.D      = 0.18;             % [M] duct/prop diameter               [m]
p.CT     = 0.10;             % [M] thrust coefficient (bench test)  [-]
p.CQ     = 0.008;            % [M] torque coefficient (bench test)  [-]
p.kappa  = (p.CQ/p.CT)*p.D;  %     reaction-torque arm  Q = kappa*T [m]
p.n_max  = 200;              % [E] max prop speed                   [rev/s]

% --- motor / inverter (rim-driven PMSM, first-order thrust model) ---
p.tau_m  = 0.03;             % [E] thrust/motor time constant       [s]
p.T_min  = 1.0;              % [E] min usable thrust per motor      [N]
p.T_max  = 14.0;             % [E] max thrust per motor             [N]
p.Kt     = 0.02;             % [M] motor torque constant            [N m/A]
p.Rphase = 0.08;             % [M] winding resistance (per phase)   [ohm]

% --- two-axis gimbal ---
p.ang_max  = deg2rad(35);    % [A] gimbal travel (both axes)        [rad]
p.rate_max = deg2rad(300);   % [E] gimbal slew-rate limit           [rad/s]
p.tau_g    = 0.03;           % [E] gimbal closed-loop time constant [s]
p.J_gimbal = 2e-4;           % [E] gimbal+prop inertia about axis   [kg m^2]

% --- rotor spin (for gyroscopic term) ---
p.J_rotor  = 1.5e-4;         % [E] rotor polar inertia              [kg m^2]
p.sig1     = +1;             %     spin handedness, propulsor 1
p.sig2     = -1;             %     spin handedness, propulsor 2 (counter-rot)

% --- battery ---
p.V_batt   = 22.2;           % [A] 6S nominal                       [V]
p.Q_batt   = 5.0;            % [A] capacity                         [Ah]

% --- hover trim (derived) ---
p.T0 = p.m*p.g/2;            % per-propulsor hover thrust           [N]
p.u_hover = [p.T0;0;0; p.T0;0;0];

% --- controller gains (per-axis, from PoC) ---
p.Kp_pos = 4.0;  p.Kd_pos = 3.5;
p.Kp_att = 40.0; p.Kd_att = 12.0;

% --- build hover effectiveness matrix B (6x6) by finite difference ---
p.B = twin_vtol_effectiveness(p);
p.Binv = pinv(p.B);
fprintf('rank(B) = %d/6, cond(B) = %.1f\n', rank(p.B), cond(p.B));

assignin('base','p',p);

function B = twin_vtol_effectiveness(p)
    u0 = p.u_hover; B = zeros(6,6); eps = 1e-6;
    W0 = twin_vtol_wrench(u0,p);
    for j = 1:6
        du = u0; du(j) = du(j)+eps;
        B(:,j) = (twin_vtol_wrench(du,p)-W0)/eps;
    end
end

function W = twin_vtol_wrench(u,p)
    T1=u(1);a1=u(2);b1=u(3);T2=u(4);a2=u(5);b2=u(6);
    n1=[sin(a1);cos(a1)*sin(b1);-cos(a1)*cos(b1)];
    n2=[sin(a2);cos(a2)*sin(b2);-cos(a2)*cos(b2)];
    F1=T1*n1; F2=T2*n2;
    Q1=-p.sig1*p.kappa*F1; Q2=-p.sig2*p.kappa*F2;
    M1=cross(p.r1,F1)+Q1; M2=cross(p.r2,F2)+Q2;
    W=[F1+F2; M1+M2];
end
