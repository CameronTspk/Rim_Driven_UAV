"""
Twin Rim-Driven Vectored-Thrust UAV + REAR TAIL FAN -- 6-DOF PoC.

Change from baseline: add a THIRD rim-driven ducted fan at the tail
(behind the CG) with a SINGLE rotational DOF about the body y-axis, to
cure the pitch-authority weakness found in the twin-only layout.

Why this should work: in the twin layout both propulsors sit on the
lateral (+/-y) baseline, so the pitch moment arm is exactly zero and
pitch could only be scraped from propeller reaction torque. A fan aft of
the CG at x = -l_t finally gives pitch a real lever arm l_t.

Frames: body axes (x fwd, y right, z DOWN), NED inertial.
Attitude: unit quaternion q=[w,x,y,z] mapping BODY -> NED.

Actuators now (8): u = [T1,a1,b1,  T2,a2,b2,  T3, dlt]
  T1,T2      main fan thrusts            a,b: main fore/aft & lateral tilt
  T3         tail fan thrust             dlt: tail tilt about +y (x-z plane)
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

np.set_printoptions(precision=4, suppress=True)

# =====================================================================
# 1. REFERENCE VEHICLE  (baseline + tail fan; masses/inertia re-estimated)
# =====================================================================
g   = 9.81
m   = 1.75                      # was 1.50; +0.25 kg tail fan assembly
d   = 0.30                      # main half-span along +/-y [m]
D   = 0.18                      # duct/prop diameter [m]
CT  = 0.10                      # thrust coefficient (placeholder, MEASURE)
CQ  = 0.008                     # torque coefficient (placeholder, MEASURE)
kappa = (CQ/CT) * D             # reaction-torque arm per unit thrust [m]

l_t = 0.35                      # tail arm: tail fan is at x = -l_t [m]
h_t = 0.00                      # tail thrust-line vertical offset from CG [m]

# ---- nominal thrust split (mains + lifting tail) --------------------
W_tot = m*g
T3_0  = 4.0                              # nominal tail thrust (up) [N]
T_m0  = (W_tot - T3_0)/2.0               # nominal per-main thrust [N]
# longitudinal trim: mains must sit forward of CG so pitch balances.
#   2*T_m0*a = l_t*T3_0   ->   a = l_t*T3_0/(2*T_m0)
a = l_t*T3_0/(2.0*T_m0)                  # main fore offset (+x) [m]

# inertia re-estimate (tail mass on centerline at x=-l_t; mains shifted +a)
m_tail, m_main = 0.25, 0.20
Ixx = 0.07                                        # tail on centerline: ~0
Iyy = 0.02 + m_tail*l_t**2 + 2*m_main*a**2        # grows the most
Izz = 0.075 + m_tail*l_t**2 + 2*m_main*a**2
I    = np.diag([Ixx, Iyy, Izz])
Iinv = np.linalg.inv(I)

# fan positions relative to CG
r1 = np.array([ a,  d, 0.0])     # right main (+y), forward
r2 = np.array([ a, -d, 0.0])     # left  main (-y), forward
r3 = np.array([-l_t, 0.0, h_t])  # tail (-x)
sig1, sig2, sig3 = +1.0, -1.0, +1.0     # spin handedness (counter-rot mains)

# actuator limits / dynamics                     (8-vector ordering as above)
T_min,  T_max  = 1.0, 14.0        # main thrust bounds [N]
Tt_min, Tt_max = 0.5, 10.0        # tail thrust bounds [N]
ang_max  = np.deg2rad(35.0)       # main gimbal travel [rad]
dlt_max  = np.deg2rad(45.0)       # tail tilt travel [rad]
rate_max = np.deg2rad(300.0)      # gimbal/tilt slew rate [rad/s]
Tdot_max = 200.0                  # thrust slew [N/s]
tau_g, tau_m = 0.03, 0.03         # 1st-order actuator lags [s]

# =====================================================================
# 2. ACTUATOR PHYSICS: thrust vectors + wrench (3 fans)
# =====================================================================
def main_dir(alpha, beta):
    """Main-fan unit thrust dir. alpha: fore/aft tilt(+x); beta: lateral(+y)."""
    ca, sa = np.cos(alpha), np.sin(alpha)
    cb, sb = np.cos(beta),  np.sin(beta)
    return np.array([sa, ca*sb, -ca*cb])

def tail_dir(dlt):
    """Tail-fan unit thrust dir: nominal up (0,0,-1), tilt dlt about +y
       swings it in the x-z plane -> (-sin dlt, 0, -cos dlt)."""
    return np.array([-np.sin(dlt), 0.0, -np.cos(dlt)])

def wrench(u):
    T1,a1,b1, T2,a2,b2, T3,dl = u
    n1, n2, n3 = main_dir(a1,b1), main_dir(a2,b2), tail_dir(dl)
    F1, F2, F3 = T1*n1, T2*n2, T3*n3
    Q1 = -sig1*kappa*F1;  Q2 = -sig2*kappa*F2;  Q3 = -sig3*kappa*F3
    M1 = np.cross(r1,F1)+Q1
    M2 = np.cross(r2,F2)+Q2
    M3 = np.cross(r3,F3)+Q3
    return np.concatenate([F1+F2+F3, M1+M2+M3])

u_hover = np.array([T_m0,0,0,  T_m0,0,0,  T3_0, 0.0])
W0 = wrench(u_hover)

# =====================================================================
# 3. EFFECTIVENESS MATRIX (6x8) + RANK / CONDITION / AUTHORITY
# =====================================================================
NU = 8
def effectiveness(u, eps=1e-6):
    B = np.zeros((6, NU)); W = wrench(u)
    for j in range(NU):
        du = u.copy(); du[j] += eps
        B[:, j] = (wrench(du) - W)/eps
    return B

B = effectiveness(u_hover)
Bpinv = B.T @ np.linalg.inv(B @ B.T)         # min-norm right pseudoinverse
S_ = np.linalg.svd(B, compute_uv=False)
rank = np.linalg.matrix_rank(B, tol=1e-6)
cond = S_[0]/S_[-1]

def max_axis_moment():
    dT_m = min(T_max - T_m0, T_m0 - T_min)
    Mx = abs(d * (2*dT_m))                          # roll: differential main thrust
    Mz = abs(d * T_m0 * 2*np.sin(ang_max))          # yaw: differential fore/aft tilt
    dT_t = min(Tt_max - T3_0, T3_0 - Tt_min)
    My = abs(l_t * dT_t)                            # PITCH: tail thrust modulation @ arm l_t
    return Mx, My, Mz

Mx_max, My_max, Mz_max = max_axis_moment()
ang_acc = np.array([Mx_max/Ixx, My_max/Iyy, Mz_max/Izz])

print("="*68)
print("TAIL-FAN CONFIG -- EFFECTIVENESS / CONTROLLABILITY (hover lin.)")
print("="*68)
print(f"mass {m:.2f} kg | tail arm l_t {l_t:.2f} m | main fore offset a {a*1000:5.1f} mm")
print(f"nominal thrust: mains {T_m0:.2f} N each, tail {T3_0:.2f} N")
print(f"inertia diag [Ixx Iyy Izz] = [{Ixx:.3f} {Iyy:.3f} {Izz:.3f}] kg m^2")
print(f"kappa (reaction arm) {kappa*1000:.2f} mm")
print(f"\nHover wrench [Fx Fy Fz Mx My Mz]: {W0}")
print("\nB = d[Fx Fy Fz Mx My Mz]/d[T1 a1 b1  T2 a2 b2  T3 dlt]  (6x8):")
print(B)
print(f"\nrank(B)          : {rank} / 6  -> {'FULL RANK (6-DOF reachable)' if rank==6 else 'RANK DEFICIENT'}")
print(f"singular values  : {S_}")
print(f"condition number : {cond:.1f}   (baseline twin-only was ~69)")
print("\nMax moment authority per axis (vertical force held = weight):")
print(f"  Roll  Mx_max = {Mx_max:6.3f} N m -> {ang_acc[0]:6.1f} rad/s^2")
print(f"  Pitch My_max = {My_max:6.3f} N m -> {ang_acc[1]:6.1f} rad/s^2   <-- was 6.1")
print(f"  Yaw   Mz_max = {Mz_max:6.3f} N m -> {ang_acc[2]:6.1f} rad/s^2")
print(f"\nPitch/roll authority ratio: {ang_acc[1]/ang_acc[0]:.2f} "
      f"(was 0.11);  pitch/yaw: {ang_acc[1]/ang_acc[2]:.2f} (was 0.18)")

# pitch effectiveness split: thrust-modulation vs tilt (columns 6,7 of B, row My=4)
print("\nPitch-moment effectiveness (row My):")
print(f"  d(My)/d(T3)  = {B[4,6]:+.4f} N m per N     (tail thrust modulation -- primary)")
print(f"  d(My)/d(dlt) = {B[4,7]:+.4f} N m per rad   (tail tilt about y; grows with h_t)")

# =====================================================================
# 4. QUATERNION + 6-DOF DYNAMICS
# =====================================================================
def qmul(a,b):
    w1,x1,y1,z1=a; w2,x2,y2,z2=b
    return np.array([w1*w2-x1*x2-y1*y2-z1*z2, w1*x2+x1*w2+y1*z2-z1*y2,
                     w1*y2-x1*z2+y1*w2+z1*x2, w1*z2+x1*y2-y1*x2+z1*w2])
def qconj(a): return np.array([a[0],-a[1],-a[2],-a[3]])
def q2R(q):
    w,x,y,z=q
    return np.array([[1-2*(y*y+z*z),2*(x*y-w*z),2*(x*z+w*y)],
                     [2*(x*y+w*z),1-2*(x*x+z*z),2*(y*z-w*x)],
                     [2*(x*z-w*y),2*(y*z+w*x),1-2*(x*x+y*y)]])
def q2euler(q):
    w,x,y,z=q
    return np.array([np.arctan2(2*(w*x+y*z),1-2*(x*x+y*y)),
                     np.arcsin(np.clip(2*(w*y-z*x),-1,1)),
                     np.arctan2(2*(w*z+x*y),1-2*(y*y+z*z))])

def rb_deriv(x, W):
    v=x[3:6]; q=x[6:10]; w=x[10:13]
    R=q2R(q)
    F_ned = R@W[0:3] + np.array([0,0,m*g])
    return np.concatenate([v, F_ned/m, 0.5*qmul(q,np.array([0,*w])),
                           Iinv@(W[3:6]-np.cross(w, I@w))])

# =====================================================================
# 5. CASCADED CONTROLLER + ALLOCATION (pseudoinverse over 8 actuators)
# =====================================================================
Kp_pos, Kd_pos = 4.0, 3.5
Kp_att, Kd_att = 40.0, 12.0
lo = np.array([T_min,-ang_max,-ang_max, T_min,-ang_max,-ang_max, Tt_min,-dlt_max])
hi = np.array([T_max, ang_max, ang_max, T_max, ang_max, ang_max, Tt_max, dlt_max])

def controller(x, ref):
    p=x[0:3]; v=x[3:6]; q=x[6:10]; w=x[10:13]
    p_des,v_des,q_des = ref
    R=q2R(q)
    a_des = Kp_pos*(p_des-p)+Kd_pos*(v_des-v)
    F_body_des = R.T @ (m*a_des - np.array([0,0,m*g]))
    qerr = qmul(qconj(q_des), q)
    if qerr[0]<0: qerr=-qerr
    ang_acc_des = -Kp_att*(2*qerr[1:4]) - Kd_att*w
    M_des = I@ang_acc_des + np.cross(w, I@w)
    W_des = np.concatenate([F_body_des, M_des])
    u = u_hover + Bpinv @ (W_des - W0)
    return np.clip(u, lo, hi)

tau = np.array([tau_m,tau_g,tau_g, tau_m,tau_g,tau_g, tau_m,tau_g])
rmx = np.array([Tdot_max,rate_max,rate_max, Tdot_max,rate_max,rate_max,
                Tdot_max,rate_max])
def actuator_step(u_act, u_cmd, dt):
    tgt = u_act + (u_cmd-u_act)*(dt/tau)
    du  = np.clip(tgt-u_act, -rmx*dt, rmx*dt)
    return np.clip(u_act+du, lo, hi)

def simulate(t_end, ref_fn, x0, dt=0.002):
    N=int(t_end/dt)
    X=np.zeros((N,13)); U=np.zeros((N,NU)); Tl=np.zeros(N)
    x=x0.copy(); u_act=u_hover.copy()
    for k in range(N):
        t=k*dt; ref=ref_fn(t)
        u_act=actuator_step(u_act, controller(x,ref), dt)
        W=wrench(u_act)
        k1=rb_deriv(x,W); k2=rb_deriv(x+0.5*dt*k1,W)
        k3=rb_deriv(x+0.5*dt*k2,W); k4=rb_deriv(x+dt*k3,W)
        x=x+(dt/6)*(k1+2*k2+2*k3+k4)
        x[6:10]/=np.linalg.norm(x[6:10])
        X[k]=x; U[k]=u_act; Tl[k]=t
    return Tl,X,U

def x_init(p=(0,0,0), v=(0,0,0), q=(1,0,0,0), w=(0,0,0)):
    return np.concatenate([p,v,q,w]).astype(float)
level=(1,0,0,0)
hold = lambda pd: (lambda t:(np.array(pd,float),np.zeros(3),np.array(level,float)))

tA,XA,UA    = simulate(4.0, hold((0,0,0)), x_init())
tB1,XB1,UB1 = simulate(3.0, hold((0,0,0)), x_init(w=(2.0,0,0)))   # roll 2
tB2,XB2,UB2 = simulate(3.0, hold((0,0,0)), x_init(w=(0,2.0,0)))   # pitch 2
tB3,XB3,UB3 = simulate(3.0, hold((0,0,0)), x_init(w=(0,4.0,0)))   # pitch 4 (the old near-tumble)
tC,XC,UC    = simulate(6.0, hold((1,0,0)), x_init())              # +1 m x step

eA =np.array([q2euler(q) for q in XA[:,6:10]])
e1 =np.array([q2euler(q) for q in XB1[:,6:10]])
e2 =np.array([q2euler(q) for q in XB2[:,6:10]])
e3 =np.array([q2euler(q) for q in XB3[:,6:10]])
eC =np.array([q2euler(q) for q in XC[:,6:10]])
def settle(t,s,tol=np.deg2rad(2.0)):
    idx=np.where(np.abs(s)>tol)[0]; return t[idx[-1]] if len(idx) else 0.0

print("\n"+"="*68); print("CLOSED-LOOP RESULTS  (tail-fan config)"); print("="*68)
print(f"Hover: max attitude {np.rad2deg(np.abs(eA).max()):.3f} deg, "
      f"drift {np.abs(XA[:,0:3]).max()*1000:.2f} mm")
print(f"Roll  disturb 2 rad/s : peak {np.rad2deg(np.abs(e1[:,0]).max()):5.1f} deg, "
      f"settle {settle(tB1,e1[:,0]):.2f} s")
print(f"Pitch disturb 2 rad/s : peak {np.rad2deg(np.abs(e2[:,1]).max()):5.1f} deg, "
      f"settle {settle(tB2,e2[:,1]):.2f} s   (twin-only was 25.9 deg / 1.02 s)")
print(f"Pitch disturb 4 rad/s : peak {np.rad2deg(np.abs(e3[:,1]).max()):5.1f} deg, "
      f"settle {settle(tB3,e3[:,1]):.2f} s   (twin-only was 89.9 deg / 2.95 s, near-tumble)")
tsat = np.mean(np.abs(UB3[:,7])>0.98*dlt_max)*100
ttsat= np.mean((UB3[:,6]>0.98*Tt_max)|(UB3[:,6]<1.02*Tt_min))*100
print(f"Tail-tilt saturation during hard pitch recovery: {tsat:.1f}% of samples")
print(f"Tail-thrust saturation during hard pitch recovery: {ttsat:.1f}% of samples")
print(f"Position step +1 m: settle {settle(tC,XC[:,0]-1.0,0.02):.2f} s, "
      f"max tilt {np.rad2deg(np.abs(eC[:,1]).max()):.2f} deg")

# =====================================================================
# 6. PLOTS
# =====================================================================
fig, ax = plt.subplots(2, 2, figsize=(12, 8))

ax[0,0].plot(tB1, np.rad2deg(e1[:,0]), label="roll, 2 rad/s")
ax[0,0].plot(tB2, np.rad2deg(e2[:,1]), label="pitch, 2 rad/s (tail fan)")
ax[0,0].plot(tB3, np.rad2deg(e3[:,1]), '--', label="pitch, 4 rad/s (tail fan)")
ax[0,0].axhline(0,color='k',lw=0.5)
ax[0,0].set_title("Disturbance recovery WITH tail fan\n(pitch now recovers like roll; cf. near-tumble before)")
ax[0,0].set_xlabel("time [s]"); ax[0,0].set_ylabel("angle [deg]"); ax[0,0].legend(); ax[0,0].grid(alpha=0.3)

names=["Roll","Pitch","Yaw"]
old=np.array([54.5,6.1,33.8]); new=ang_acc
xpos=np.arange(3); wd=0.38
ax[0,1].bar(xpos-wd/2, old, wd, label="twin only", color="#bbb")
ax[0,1].bar(xpos+wd/2, new, wd, label="with tail fan", color=["#2a7","#d33","#27a"])
ax[0,1].set_xticks(xpos); ax[0,1].set_xticklabels(names)
ax[0,1].set_title("Angular-acceleration authority per axis\ntwin-only vs tail-fan")
ax[0,1].set_ylabel("max ang. accel [rad/s$^2$]"); ax[0,1].legend(); ax[0,1].grid(alpha=0.3,axis='y')
for xp,v in zip(xpos+wd/2,new): ax[0,1].text(xp,v,f"{v:.0f}",ha='center',va='bottom')
for xp,v in zip(xpos-wd/2,old): ax[0,1].text(xp,v,f"{v:.0f}",ha='center',va='bottom',color='#666')

ax[1,0].plot(tB3, UB3[:,6], label=r"tail thrust $T_3$ [N]")
ax[1,0].axhline(Tt_max,color='r',ls=':',label="thrust limit")
ax[1,0].axhline(Tt_min,color='r',ls=':')
ax[1,0].plot(tB3, np.rad2deg(UB3[:,7])/10.0, label=r"tail tilt $\delta$ [deg/10]")
ax[1,0].set_title("Tail actuator during hard (4 rad/s) pitch recovery\n(primary pitch lever = tail thrust)")
ax[1,0].set_xlabel("time [s]"); ax[1,0].set_ylabel("thrust [N] / tilt"); ax[1,0].legend(); ax[1,0].grid(alpha=0.3)

ax[1,1].plot(tC, XC[:,0], label="x position")
ax[1,1].plot(tC, np.rad2deg(eC[:,1]), label="pitch angle")
ax[1,1].axhline(1.0,color='k',ls=':',lw=0.8)
ax[1,1].set_title("Position step +1 m in x\n(still translates at ~0 tilt: fully-actuated advantage kept)")
ax[1,1].set_xlabel("time [s]"); ax[1,1].set_ylabel("x [m] / pitch [deg]"); ax[1,1].legend(); ax[1,1].grid(alpha=0.3)

plt.tight_layout()
plt.savefig("/home/claude/tail_poc_results.png", dpi=120, bbox_inches="tight")
print("\nSaved: tail_poc_results.png")
