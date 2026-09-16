"""
Twin Rim-Driven Vectored-Thrust UAV -- Proof-of-Concept 6-DOF simulation.

Purpose (per the brief):
  1. Build a clearly-defined reference vehicle.
  2. Derive the control-effectiveness (wrench) matrix and test its RANK and
     CONDITION NUMBER -- i.e. do 6 actuators really give 6 controllable DOF?
  3. Run closed-loop tests (hover, roll vs pitch disturbance, position step)
     to confirm the plant is controllable BEFORE committing to Simulink/PX4.

Frames: aerospace body axes (x fwd, y right, z DOWN) and NED inertial.
Attitude: unit quaternion q = [w,x,y,z] mapping BODY -> NED.
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

np.set_printoptions(precision=4, suppress=True)

# =====================================================================
# 1. REFERENCE VEHICLE PARAMETERS  (realistic small-UAV placeholders)
# =====================================================================
g      = 9.81
m      = 1.50                      # vehicle mass [kg]
d      = 0.30                      # half-span: propulsor offset along +/-y [m]
D      = 0.18                      # duct/prop diameter [m]
CT     = 0.10                      # thrust coefficient (placeholder)
CQ     = 0.008                     # torque coefficient (placeholder)
kappa  = (CQ / CT) * D             # reaction torque per unit thrust [m]
I      = np.diag([0.07, 0.02, 0.075])   # inertia Ixx,Iyy,Izz [kg m^2]
Iinv   = np.linalg.inv(I)
T0     = m * g / 2.0               # per-propulsor hover thrust [N]

# propulsor positions relative to CG (lateral baseline)
r1 = np.array([0.0,  d, 0.0])      # right propulsor (+y)
r2 = np.array([0.0, -d, 0.0])      # left  propulsor (-y)
sig1, sig2 = +1.0, -1.0            # counter-rotating spin handedness

# actuator limits / dynamics
T_min, T_max = 1.0, 14.0           # per-motor thrust bounds [N]  (T/W_max ~ 1.9)
ang_max  = np.deg2rad(35.0)        # gimbal travel [rad]
rate_max = np.deg2rad(300.0)       # gimbal slew rate [rad/s]
Tdot_max = 200.0                   # thrust slew [N/s]
tau_g    = 0.03                    # gimbal 1st-order lag [s]
tau_m    = 0.03                    # motor thrust lag [s]

# =====================================================================
# 2. ACTUATOR PHYSICS: thrust vector + wrench
# =====================================================================
def thrust_dir(alpha, beta):
    """Unit thrust direction in body axes. alpha: fore/aft tilt (+x),
       beta: lateral tilt (+y). Nominal (0,0) -> (0,0,-1) = up."""
    ca, sa = np.cos(alpha), np.sin(alpha)
    cb, sb = np.cos(beta),  np.sin(beta)
    return np.array([sa, ca * sb, -ca * cb])

def wrench(u):
    """u = [T1, a1, b1, T2, a2, b2] -> W = [Fx,Fy,Fz, Mx,My,Mz] in body axes."""
    T1, a1, b1, T2, a2, b2 = u
    n1 = thrust_dir(a1, b1); n2 = thrust_dir(a2, b2)
    F1 = T1 * n1;            F2 = T2 * n2
    Q1 = -sig1 * kappa * F1; Q2 = -sig2 * kappa * F2      # reaction torque
    M1 = np.cross(r1, F1) + Q1
    M2 = np.cross(r2, F2) + Q2
    return np.concatenate([F1 + F2, M1 + M2])

u_hover = np.array([T0, 0.0, 0.0, T0, 0.0, 0.0])
W0 = wrench(u_hover)

# =====================================================================
# 3. CONTROL-EFFECTIVENESS MATRIX + RANK / CONDITION / AUTHORITY
# =====================================================================
def effectiveness(u, eps=1e-6):
    B = np.zeros((6, 6)); W = wrench(u)
    for j in range(6):
        du = u.copy(); du[j] += eps
        B[:, j] = (wrench(du) - W) / eps
    return B

B    = effectiveness(u_hover)
Binv = np.linalg.inv(B)
U_, S_, Vt_ = np.linalg.svd(B)
rank = np.linalg.matrix_rank(B, tol=1e-6)
cond = S_[0] / S_[-1]

def max_axis_moment():
    """Representative max moment about each body axis while keeping total
       vertical thrust = weight (so the vehicle doesn't fall)."""
    # Roll: differential thrust, sum held at 2*T0
    dT   = min(T_max - T0, T0 - T_min)
    Mx   = abs(d * (2 * dT))
    # Yaw: differential fore/aft tilt at hover thrust
    Mz   = abs(d * T0 * 2 * np.sin(ang_max))
    # Pitch: ONLY differential lateral tilt via reaction torque
    My   = abs(kappa * T0 * 2 * np.sin(ang_max))
    return Mx, My, Mz

Mx_max, My_max, Mz_max = max_axis_moment()
ang_acc = np.array([Mx_max / I[0, 0], My_max / I[1, 1], Mz_max / I[2, 2]])

print("=" * 66)
print("EFFECTIVENESS / CONTROLLABILITY ANALYSIS  (hover linearization)")
print("=" * 66)
print(f"kappa (reaction arm)         : {kappa*1000:.2f} mm")
print(f"Hover wrench [Fx Fy Fz Mx My Mz]: {W0}")
print("\nB = d[Fx Fy Fz Mx My Mz]/d[T1 a1 b1 T2 a2 b2]:")
print(B)
print(f"\nrank(B)            : {rank} / 6   -> {'FULL RANK (6-DOF reachable)' if rank==6 else 'RANK DEFICIENT'}")
print(f"singular values    : {S_}")
print(f"condition number   : {cond:.1f}")
print("\nMax moment authority about each axis (vertical force held = weight):")
print(f"  Roll  Mx_max = {Mx_max:6.3f} N m -> {ang_acc[0]:6.1f} rad/s^2")
print(f"  Pitch My_max = {My_max:6.3f} N m -> {ang_acc[1]:6.1f} rad/s^2   <-- WEAK AXIS")
print(f"  Yaw   Mz_max = {Mz_max:6.3f} N m -> {ang_acc[2]:6.1f} rad/s^2")
print(f"\nPitch authority is {ang_acc[0]/ang_acc[1]:.1f}x weaker than roll, "
      f"{ang_acc[2]/ang_acc[1]:.1f}x weaker than yaw.")

# =====================================================================
# 4. QUATERNION HELPERS + 6-DOF DYNAMICS
# =====================================================================
def qmul(a, b):
    w1,x1,y1,z1 = a; w2,x2,y2,z2 = b
    return np.array([w1*w2-x1*x2-y1*y2-z1*z2,
                     w1*x2+x1*w2+y1*z2-z1*y2,
                     w1*y2-x1*z2+y1*w2+z1*x2,
                     w1*z2+x1*y2-y1*x2+z1*w2])
def qconj(a): return np.array([a[0], -a[1], -a[2], -a[3]])
def q2R(q):
    w,x,y,z = q
    return np.array([
        [1-2*(y*y+z*z), 2*(x*y-w*z),   2*(x*z+w*y)],
        [2*(x*y+w*z),   1-2*(x*x+z*z), 2*(y*z-w*x)],
        [2*(x*z-w*y),   2*(y*z+w*x),   1-2*(x*x+y*y)]])
def q2euler(q):
    w,x,y,z = q
    roll  = np.arctan2(2*(w*x+y*z), 1-2*(x*x+y*y))
    pitch = np.arcsin(np.clip(2*(w*y-z*x), -1, 1))
    yaw   = np.arctan2(2*(w*z+x*y), 1-2*(y*y+z*z))
    return np.array([roll, pitch, yaw])

def rb_deriv(x, W):
    p = x[0:3]; v = x[3:6]; q = x[6:10]; w = x[10:13]
    F_body = W[0:3]; M_body = W[3:6]
    R = q2R(q)
    F_ned = R @ F_body + np.array([0, 0, m*g])       # + gravity (NED, down = +z)
    pdot = v
    vdot = F_ned / m
    qdot = 0.5 * qmul(q, np.array([0, *w]))
    wdot = Iinv @ (M_body - np.cross(w, I @ w))
    return np.concatenate([pdot, vdot, qdot, wdot])

# =====================================================================
# 5. CASCADED CONTROLLER + CONTROL ALLOCATION
# =====================================================================
Kp_pos, Kd_pos = 4.0, 3.5
Kp_att, Kd_att = 40.0, 12.0

def controller(x, ref):
    p = x[0:3]; v = x[3:6]; q = x[6:10]; w = x[10:13]
    p_des, v_des, q_des = ref
    R = q2R(q)
    # outer position loop -> desired body force (fully-actuated: no tilt needed)
    a_des = Kp_pos*(p_des - p) + Kd_pos*(v_des - v)
    F_thr_ned = m*a_des - np.array([0, 0, m*g])      # cancel gravity, add accel
    F_body_des = R.T @ F_thr_ned
    # attitude loop -> desired moment
    qerr = qmul(qconj(q_des), q)
    if qerr[0] < 0: qerr = -qerr
    e_ang = 2.0 * qerr[1:4]
    ang_acc_des = -Kp_att*e_ang - Kd_att*w
    M_des = I @ ang_acc_des + np.cross(w, I @ w)
    # allocation via matrix inverse of effectiveness matrix
    W_des = np.concatenate([F_body_des, M_des])
    u_cmd = u_hover + Binv @ (W_des - W0)
    u_cmd[[0, 3]]       = np.clip(u_cmd[[0, 3]], T_min, T_max)
    u_cmd[[1, 2, 4, 5]] = np.clip(u_cmd[[1, 2, 4, 5]], -ang_max, ang_max)
    return u_cmd

def actuator_step(u_act, u_cmd, dt):
    """First-order lag + rate + position saturation on each actuator."""
    tgt = u_act + (u_cmd - u_act) * (dt / np.array([tau_m, tau_g, tau_g,
                                                    tau_m, tau_g, tau_g]))
    du  = tgt - u_act
    rmax = np.array([Tdot_max, rate_max, rate_max,
                     Tdot_max, rate_max, rate_max]) * dt
    du   = np.clip(du, -rmax, rmax)
    u    = u_act + du
    u[[0, 3]]       = np.clip(u[[0, 3]], T_min, T_max)
    u[[1, 2, 4, 5]] = np.clip(u[[1, 2, 4, 5]], -ang_max, ang_max)
    return u

def simulate(t_end, ref_fn, x0, dt=0.002):
    N = int(t_end/dt)
    X = np.zeros((N, 13)); U = np.zeros((N, 6)); Tlog = np.zeros(N)
    x = x0.copy(); u_act = u_hover.copy()
    for k in range(N):
        t = k*dt
        ref = ref_fn(t)
        u_cmd = controller(x, ref)
        u_act = actuator_step(u_act, u_cmd, dt)
        W = wrench(u_act)
        # RK4 on rigid body with fixed wrench over dt
        k1 = rb_deriv(x, W)
        k2 = rb_deriv(x + 0.5*dt*k1, W)
        k3 = rb_deriv(x + 0.5*dt*k2, W)
        k4 = rb_deriv(x + dt*k3, W)
        x  = x + (dt/6.0)*(k1 + 2*k2 + 2*k3 + k4)
        x[6:10] /= np.linalg.norm(x[6:10])            # renormalize quaternion
        X[k] = x; U[k] = u_act; Tlog[k] = t
    return Tlog, X, U

def x_init(p=(0,0,0), v=(0,0,0), q=(1,0,0,0), w=(0,0,0)):
    return np.concatenate([p, v, q, w]).astype(float)

level = (1,0,0,0)
hold  = lambda p_des: (lambda t: (np.array(p_des,float), np.zeros(3), np.array(level,float)))

# ---- Test A: hover trim hold ----
tA, XA, UA = simulate(4.0, hold((0,0,0)), x_init())

# ---- Test B: roll vs pitch disturbance recovery (initial body rate 2 rad/s) ----
tB1, XB1, UB1 = simulate(3.0, hold((0,0,0)), x_init(w=(2.0,0,0)))   # roll
tB2, XB2, UB2 = simulate(3.0, hold((0,0,0)), x_init(w=(0,2.0,0)))   # pitch
tB3, XB3, UB3 = simulate(3.0, hold((0,0,0)), x_init(w=(0,4.0,0)))   # pitch (harder)

# ---- Test C: position step +1 m in x (translate at level attitude) ----
tC, XC, UC = simulate(6.0, hold((1,0,0)), x_init())

eulB1 = np.array([q2euler(q) for q in XB1[:, 6:10]])
eulB2 = np.array([q2euler(q) for q in XB2[:, 6:10]])
eulB3 = np.array([q2euler(q) for q in XB3[:, 6:10]])
eulA  = np.array([q2euler(q) for q in XA[:,  6:10]])

def settle_time(t, sig, tol=np.deg2rad(2.0)):
    idx = np.where(np.abs(sig) > tol)[0]
    return t[idx[-1]] if len(idx) else 0.0

print("\n" + "="*66)
print("CLOSED-LOOP RESULTS")
print("="*66)
print(f"Hover: max attitude excursion = {np.rad2deg(np.abs(eulA).max()):.3f} deg, "
      f"max drift = {np.abs(XA[:,0:3]).max()*1000:.2f} mm  (trim OK)")
print(f"Roll  disturb 2 rad/s : peak {np.rad2deg(np.abs(eulB1[:,0]).max()):5.1f} deg, "
      f"settle {settle_time(tB1, eulB1[:,0]):.2f} s")
print(f"Pitch disturb 2 rad/s : peak {np.rad2deg(np.abs(eulB2[:,1]).max()):5.1f} deg, "
      f"settle {settle_time(tB2, eulB2[:,1]):.2f} s")
print(f"Pitch disturb 4 rad/s : peak {np.rad2deg(np.abs(eulB3[:,1]).max()):5.1f} deg, "
      f"settle {settle_time(tB3, eulB3[:,1]):.2f} s")
# gimbal saturation fraction during pitch recovery
bsat = np.mean(np.abs(UB3[:, [2,5]]) > 0.98*ang_max) * 100
print(f"Lateral-gimbal saturation during hard pitch recovery: {bsat:.1f}% of samples")
print(f"Position step +1 m: settle {settle_time(tC, XC[:,0]-1.0, 0.02):.2f} s, "
      f"max tilt during move = {np.rad2deg(np.abs(np.array([q2euler(q) for q in XC[:,6:10]])[:,1]).max()):.2f} deg")

# =====================================================================
# 6. PLOTS
# =====================================================================
fig, ax = plt.subplots(2, 2, figsize=(12, 8))

ax[0,0].plot(tB1, np.rad2deg(eulB1[:,0]), label="roll axis (strong)")
ax[0,0].plot(tB2, np.rad2deg(eulB2[:,1]), label="pitch axis (weak)")
ax[0,0].plot(tB3, np.rad2deg(eulB3[:,1]), '--', label="pitch axis, 2x disturbance")
ax[0,0].axhline(0, color='k', lw=0.5)
ax[0,0].set_title("Disturbance recovery: roll vs pitch\n(same 2 rad/s initial rate)")
ax[0,0].set_xlabel("time [s]"); ax[0,0].set_ylabel("angle [deg]"); ax[0,0].legend(); ax[0,0].grid(alpha=0.3)

axes_names = ["Roll","Pitch","Yaw"]
bars = ax[0,1].bar(axes_names, ang_acc, color=["#2a7","#d33","#27a"])
ax[0,1].set_title("Angular-acceleration authority per axis\n(vertical thrust held = weight)")
ax[0,1].set_ylabel("max ang. accel [rad/s$^2$]"); ax[0,1].grid(alpha=0.3, axis='y')
for b,val in zip(bars, ang_acc):
    ax[0,1].text(b.get_x()+b.get_width()/2, val, f"{val:.0f}", ha='center', va='bottom')

ax[1,0].plot(tB3, np.rad2deg(UB3[:,2]), label=r"$\beta_1$ (right)")
ax[1,0].plot(tB3, np.rad2deg(UB3[:,5]), label=r"$\beta_2$ (left)")
ax[1,0].axhline( np.rad2deg(ang_max), color='r', ls=':', label="limit")
ax[1,0].axhline(-np.rad2deg(ang_max), color='r', ls=':')
ax[1,0].set_title("Lateral-gimbal command during hard pitch recovery\n(pitch uses the SAME actuator as lateral force)")
ax[1,0].set_xlabel("time [s]"); ax[1,0].set_ylabel("gimbal $\\beta$ [deg]"); ax[1,0].legend(); ax[1,0].grid(alpha=0.3)

eulC = np.array([q2euler(q) for q in XC[:, 6:10]])
ax[1,1].plot(tC, XC[:,0], label="x position")
ax[1,1].plot(tC, np.rad2deg(eulC[:,1]), label="pitch angle")
ax[1,1].axhline(1.0, color='k', ls=':', lw=0.8)
ax[1,1].set_title("Position step +1 m in x\n(vectored: translates with ~0 tilt)")
ax[1,1].set_xlabel("time [s]"); ax[1,1].set_ylabel("x [m]  /  pitch [deg]"); ax[1,1].legend(); ax[1,1].grid(alpha=0.3)

plt.tight_layout()
plt.savefig("/home/claude/poc_results.png", dpi=120, bbox_inches="tight")
print("\nSaved: poc_results.png")
