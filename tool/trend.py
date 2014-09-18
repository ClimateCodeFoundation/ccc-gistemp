#!/usr/bin/env python

def lm1(data):
    """
    Regresses the linear model y = a + b*x and returns (a, b,
    r2) for the best fit (where r2 is the R-squared value).

    The input, *data*, should be a series of (x,y) pairs. y can
    be None in which case the pair is ignored.
    """
    sxx = sxy = syy = sx = sy = n = 0
    for (x,y) in data:
        if y is not None:
            sxx += x*x
            syy += y*y
            sx += x
            sy += y
            sxy += x * y
            n += 1
    if n < 2:
        return None,None,None
    # Make n a float. This contaminates all the subsequent divisions,
    # making them floating point divisions with floating point answers,
    # which is what we want.
    n = float(n)
    xbar = sx / n
    ybar = sy / n
    ssxx = sxx - (sx * sx) / n
    ssyy = syy - (sy * sy) / n
    ssxy = sxy - (sx * sy) / n
    if ssxx == 0:
        return None, None, None
    b = ssxy / ssxx
    a = ybar - b * xbar
    r2 = (ssxy * ssxy) / (ssxx * ssyy)
    return (a,b, r2)
