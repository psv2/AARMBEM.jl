"""

    E1func(x)

Compute the exponential integral function ``E_{1}(x)`` for a real number `x`. This function
returns zero if `x` is larger than 700, or `Inf + pi*1im` if `x` is smaller than -700. (This
definition also defines the branch taken in the complex plane for ``x < 0``.) Note that even
though ``E_{1}(x)`` is real for ``x > 0``, the result of `E1func(x)` is cast as a complex
type, based on the type of `x`, for consistency.

This is a reimplementation of code at
`https://github.com/mitmath/18S096/blob/iap2017/pset3/pset3-solutions.ipynb` by Steven G.
Johnson, based on a combination of Taylor series and continued fraction expansions.

See also: [`E1taylor`](@ref), [`E1cf`](@ref)

"""
function E1func(x::AbstractFloat)
  
    if (x < -700*one(x))
        return Inf + pi*1im;
    elseif (x > 700*one(x))
        return zero(complex(x));
    else
      
        currfunc = E1taylor;
        if (x < -30*one(x) || x > 3*one(x))
            currfunc = E1cf;
        end
        mymaxorder = round(Int, 2 + 23/(1 + ((log10(abs(x)) - 0.4)/1.0)^2)); # useful empirical rule of thumb for largest order needed for positive x
        if (x < zero(x))
            mymaxorder = round(Int, 1 + 70/(1 + (log10(abs(x)) - 1.5)^2)); # useful empirical rule of thumb for largest order needed for negative x
        end
        myE1 = currfunc(x, mymaxorder);
        if (x < zero(x))
            return (complex(myE1) + pi*1im);
        end
        return complex(myE1);
        
    end
    
end

E1func(x::Real) = E1func(float(x));

"""

    E1taylor(x, maxorder)

Compute the exponential integral function ``E_{1}(x)`` for a real number `x`, using the
Taylor series expansion up to a maximum order `maxorder`. This function works best for
``-30 \\leq x \\leq 3``.

This is a reimplementation of code at
`https://github.com/mitmath/18S096/blob/iap2017/pset3/pset3-solutions.ipynb` by Steven G.
Johnson.

See also: [`E1func`](@ref), [`E1cf`](@ref)

"""
function E1taylor(x::AbstractFloat, maxorder::Integer)

    E1 = -1*MathConstants.eulergamma - log(abs(x));
    currterm = -1*x;
    E1 = E1 - currterm;
  
    for nn=2:maxorder
      
        currterm = -1*x*(nn-1)/(nn^2) * currterm;
        E1 = E1 - currterm;
        
    end
    
    return E1;
  
end

E1taylor(x::Real, maxorder::Integer) = E1taylor(float(x), maxorder);

"""

    E1cf(x, maxorder)

Compute the exponential integral function ``E_{1}(x)`` for a real number `x`, using the
continued fraction expansion up to a maximum order `maxorder`. This function works best for
``x < -30`` or ``x > 3``.

This is a reimplementation of code at
`https://github.com/mitmath/18S096/blob/iap2017/pset3/pset3-solutions.ipynb` by Steven G.
Johnson.

See also: [`E1func`](@ref), [`E1taylor`](@ref)

"""
function E1cf(x::AbstractFloat, maxorder::Integer)
  
    currcf = x;
    for nn=maxorder:-1:1
        currcf = x + (1+nn)/currcf;
        currcf = 1 + nn/currcf;
    end
    return (exp(-1*x)/(x + 1/currcf));
    
end

E1cf(x::Real, maxorder::Integer) = E1cf(float(x), maxorder);
