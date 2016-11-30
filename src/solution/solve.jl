abstract SolutionMethod
immutable VFI <: SolutionMethod
    tol
    maximum_iter
end
VFI() = VFI(1e-10, 50_000)

immutable Newton <: SolutionMethod
    tol
    maximum_iter
    d
end
Newton() = Newton(1e-10,10, 0.5)

immutable Policy <: SolutionMethod
    tol
    maximum_iter
end
Policy() = Policy(1e-10, 10)

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
end
Poly() = Poly(1e-10, 10, 1e-10, 5, 1e-10, 10, 1e-10, 10, 0.5)
# should have
# immutable Poly
# vfi
# policy
# end
function solve(U, S, method=Poly(); verbose = false)
    V = IntegratedValueFunction(S.nX, length(U.U))
    ret = solve!(U, S, V, method; verbose = verbose)
    return V, ret
end
solve(U, S, M::SolutionMethod; verbose = false) = solve!(U, S, IntegratedValueFunction(S.nX, length(U.U)), M; verbose = verbose)

function solve!(U, S, V, method::VFI; verbose = false)
    for i = 1:method.maximum_iter
        _supnorm = supnorm(V)
        bellman!(U.U, S.F, U.β, V)
        verbose && @printf("VFI i: %i   reltol: %0.5e    tol1: %0.5e\n", i, supnorm(V)/_supnorm, supnorm(V))
        if supnorm(V) < method.tol
            P!(U, S, V)
            return i
        elseif i == method.maximum_iter
            warn("Maximum number of iterations reached without convergence.")
            return i
        end
    end
end
solve!(U, S, V; verbose = false) = solve!(U, S, V, VFI(); verbose = verbose)

function solve!(U, S, V, method::Newton; tol = 1e-10, verbose = false)
    for i = 1:method.maximum_iter
        _supnorm = supnorm(V)
        newton!(U, S, V, method.d)
        verbose && @printf("Newton i: %i   reltol: %0.5e    tol1: %0.5e\n", i, supnorm(V)/_supnorm, supnorm(V))
        supnorm(V) < method.tol && return i
        i == method.maximum_iter && return i
    end
    i
end

# This is currently not "correct".
function solve!(U, S, V, method::Poly; verbose = false)
    vfi = VFI(method.tol_vfi, method.maximum_vfi)
    newton = Newton(method.tol_newton, method.maximum_newton, method.d)
    i_vfi, i_newton, i_poly = 0, 0, 0
    local _i_vfi
    for i = 1:method.maximum_poly
        for _i_vfi = 1:vfi.maximum_iter
            _supnorm = supnorm(V)
            bellman!(U.U, S.F, U.β, V)
            verbose && @printf("VFI i: %i   reltol: %0.5e    tol1: %0.5e\n", _i_vfi, supnorm(V)/_supnorm, supnorm(V))
            if supnorm(V) < vfi.tol
                break
            end
        end
        i_vfi += _i_vfi
        supnorm(V) < 1e-12 && return [i_vfi, i_newton]#, i_pol]
        j = solve!(U, S, V, newton; verbose = verbose); i_newton += j; i_vfi += j
        i_poly += 1
        supnorm(V) < 1e-12 && return [i_vfi, i_newton]#, i_poly]
        i_vfi += solve!(U, S, V, VFI(method.tol_vfi, 1000); verbose = verbose)
        supnorm(V) < 1e-12 && return [i_vfi, i_newton]#, i_poly]
    end
    supnorm(V) < 1e-12 && return [i_vfi, i_newton]#, i_poly]
end
