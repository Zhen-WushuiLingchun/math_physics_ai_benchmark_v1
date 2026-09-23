import numpy as np, time
from td_solve import *

for dx, dt in [(0.05, 0.025), (0.025, 0.0125)]:
    t0 = time.time()
    out = run_td(0.0, 12.0, dx=dx, dt=dt, Tmax=140.0)
    print(f"dx={dx} dt={dt}: {time.time()-t0:.1f}s  C={out['C']:.6f}")
    tau, h = out['tau'], out['h']
    # signal overview
    print("  |h| at tau=9,10,11,12,15,20,30,50,80,120:",
          " ".join(f"{abs(np.interp(t,tau,h)):.3e}" for t in [9,10,11,12,15,20,30,50,80,120]))
    for w in [(20, 60), (30, 80), (40, 100), (60, 120)]:
        om, sd, n = fit_frequency(tau, h, w[0], w[1])
        print(f"  fit [{w[0]},{w[1]}]: omega = {om.real:.6f} {om.imag:+.6f}i  (sd={sd:.2e}, n={n})")
    print()
