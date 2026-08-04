#!/usr/bin/env python3
"""Is the sheared restart solenoidal w.r.t. static kx or w.r.t. kxt = kx + q*tsc*ky?"""
import sys, numpy as np
sys.path.insert(0, ".")
import compare as C

N = (128, 128, 128); LX = LY = LZ = 1.0
q = 1.5
for d, code in ((sys.argv[1], sys.argv[2]),):
    cube = C.full_cube(C.load_gallope(d + "/restart", *N) if code == "gallope"
                       else C.load_calliope(d + "/restart", *N), code, *N)
    tt, tsc = C.read_restart_time(d + "/restart")[:2]
    kx, ky, kz = C.kgrids(*N, LX, LY, LZ)
    print(f"{d}  tt={tt}  tsc={tsc}")
    for kxuse, lbl in ((kx, "static kx"),
                       (kx + q * tsc * ky, "kxt = kx + q*tsc*ky"),
                       (kx - q * tsc * ky, "kxt = kx - q*tsc*ky")):
        du = C.div_measure(cube, "u", kxuse, ky, kz)
        db = C.div_measure(cube, "b", kxuse, ky, kz)
        print(f"   {lbl:24s}  <|k.u|>={du:.3e}  <|k.b|>={db:.3e}")
