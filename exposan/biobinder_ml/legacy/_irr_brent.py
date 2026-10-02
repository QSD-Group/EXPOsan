# -*- coding: utf-8 -*-
"""
Created on Mon Feb 16 16:16:45 2026

@author: aliah
"""

import numpy as np


__all__ = (
    "solve_IRR_brent",
    )

def solve_IRR_brent(tea, IRR_min=-0.9, IRR_max=5.0, maxiter=200, xtol=1e-8):
    CF = tea.cashflow_array
    dur = tea._get_duration_array()

    if CF is None or (not np.isfinite(CF).all()):
        return np.nan, "CF_nonfinite"

    def f(r):
        if r <= -0.999999:
            return np.nan
        try:
            return float((CF / (1.0 + r) ** dur).sum())
        except Exception:
            return np.nan

    a, b = IRR_min, IRR_max
    fa, fb = f(a), f(b)

    if not np.isfinite(fa) or not np.isfinite(fb):
        return np.nan, "endpoint_nonfinite"

    if np.sign(fa) == np.sign(fb):
        return np.nan, "no_bracket"

    c, fc = a, fa
    d = e = b - a

    for _ in range(maxiter):
        if fb == 0:
            return float(b), "ok"

        if np.sign(fa) == np.sign(fb):
            a, fa = c, fc
            d = e = b - a

        if abs(fa) < abs(fb):
            c, b, a = b, a, b
            fc, fb, fa = fb, fa, fb

        tol = 2 * np.finfo(float).eps * abs(b) + xtol
        m = 0.5 * (a - b)

        if abs(m) <= tol:
            return float(b), "ok"

        if abs(e) >= tol and abs(fc) > abs(fb):
            s = fb / fc
            if a == c:
                p = 2 * m * s
                q = 1 - s
            else:
                q = fc / fa
                r = fb / fa
                p = s * (2 * m * q * (q - r) - (b - c) * (r - 1))
                q = (q - 1) * (r - 1) * (s - 1)

            if p > 0:
                q = -q
            p = abs(p)

            if 2 * p < min(3 * m * q - abs(tol * q), abs(e * q)):
                e = d
                d = p / q
            else:
                d = e = m
        else:
            d = e = m

        c, fc = b, fb
        if abs(d) > tol:
            b += d
        else:
            b += tol if m > 0 else -tol

        fb = f(b)
        if not np.isfinite(fb):
            return np.nan, "solver_fail_nonfinite"

    return np.nan, "solver_fail_maxiter"