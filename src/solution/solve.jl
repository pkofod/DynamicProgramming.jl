abstract type SolutionMethod end
immutable VFI <: SolutionMethod
    tol
    maxiter::Int
    verbose::Bool
end
VFI(;maxiter=50_000) = VFI(1e-10, maxiter, false)
VFI(tol, maxiter) = VFI(tol, maxiter, false)

immutable Newton <: SolutionMethod
    tol
    maxiter::Int
    d
    verbose::Bool
end
Newton() = Newton(1e-10,10, 0.5, false)
Newton(tol, maxiter, d) = Newton(tol, maxiter, d, false)
immutable Policy <: SolutionMethod
    tol
    maxiter::Int
    verbose::Bool
end
Policy() = Policy(1e-10, 10, false)
Policy(tol, maxiter) = Policy(tol, maxiter, false)
immutable Poly <: SolutionMethod
    tol_vfi
    maximum_vfi
    tol_policy
    maximum_policy
    tol_newton
    maximum_newton
    tol_poly
    maximum_poly
    d
    verbose::Bool
end
Poly() = Poly(1e-10, 10, 1e-10, 5, 1e-10, 10, 1e-10, 10, 0.5, false)
Poly(vfitol, maxvfi, policytol, maxpolicy, newtontol, maxnewton,polytol, maxpoly, d) =
Poly(vfitol, maxvfi, policytol, maxpolicy, newtontol, maxnewton,polytol, maxpoly, d, false)
# should have
# immutable Poly
# vfi
# policy
# end
function solve!(U, S, M::SolutionMethod=Poly())
    V = IntegratedValueFunction(S.nX, length(U.U))
    ret = solve!(U, S, V, M)
    return V, ret
end

solve!(U, S, V::AbstractValueFunction) = solve!(U, S, V, VFI())
function solve!(U, S, V::AbstractValueFunction, method::VFI)
    for i = 1:method.maxiter
        _supnorm = supnorm(V)
        bellman!(U.U, S.F, U.β, V)
        method.verbose && @printf("VFI i: %i   reltol: %0.5e    tol1: %0.5e\n", i, supnorm(V)/_supnorm, supnorm(V))
        if supnorm(V) < method.tol
            P!(U, S, V)
            return i
        elseif i == method.maxiter
            warn("Maximum number of iterations reached without convergence.")
            return i
        end
    end
    i
end

function solve!(U, S, V, method::Newton)
    for i = 1:method.maxiter
        _supnorm = supnorm(V)
        newton!(U, S, V, method.d)
        method.verbose && @printf("Newton i: %i   reltol: %0.5e    tol1: %0.5e\n", i, supnorm(V)/_supnorm, supnorm(V))
        is_close(previous(V), current(V), method.tol, method.tol) && return i
        # supnorm(V) < method.tol && return i
        i == method.maxiter && return i
    end
    i
end

function solve!(U, S, V, method::Policy)
    local i::Int
    for i = 1:method.maxiter
        _supnorm = supnorm(V)
        policy!(U, S, V)
        method.verbose && @printf("Policy i: %i   reltol: %0.5e    tol1: %0.5e\n", i, supnorm(V)/_supnorm, supnorm(V))
        if is_close(previous(V), current(V), method.tol, method.tol)
            break
        end
        if i == method.maxiter
            break
        end
    end
    return i
end

# This is currently not "correct".
function solve!(U, S, V, method::Poly)
    vfi = VFI(method.tol_vfi, method.maximum_vfi)
    newton = Newton(method.tol_newton, method.maximum_newton, method.d)
    i_vfi, i_newton, i_poly = 0, 0, 0
    local _i_vfi
    for i = 1:method.maximum_poly
        for _i_vfi = 1:vfi.maxiter
            _supnorm = supnorm(V)
            bellman!(U.U, S.F, U.β, V)
            method.verbose && @printf("VFI i: %i   reltol: %0.5e    tol1: %0.5e\n", _i_vfi, supnorm(V)/_supnorm, supnorm(V))
            if supnorm(V) < vfi.tol
                break
            end
        end
        i_vfi += _i_vfi
        supnorm(V) < 1e-12 && return [i_vfi, i_newton]#, i_pol]
        j = solve!(U, S, V, newton); i_newton += j; i_vfi += j
        i_poly += 1
        supnorm(V) < 1e-12 && return [i_vfi, i_newton]#, i_poly]
        i_vfi += solve!(U, S, V, VFI(method.tol_vfi, 1000))
        supnorm(V) < 1e-12 && return [i_vfi, i_newton]#, i_poly]
    end
    supnorm(V) < 1e-12 && return [i_vfi, i_newton]#, i_poly]
end
