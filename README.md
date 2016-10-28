# MDPTools

[![Build Status](https://travis-ci.org/pkofod/MDPTools.jl.svg?branch=master)](https://travis-ci.org/pkofod/MDPTools.jl)

[![Coverage Status](https://coveralls.io/repos/pkofod/MDPTools.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/pkofod/MDPTools.jl?branch=master)

[![codecov.io](http://codecov.io/github/pkofod/MDPTools.jl/coverage.svg?branch=master)](http://codecov.io/github/pkofod/MDPTools.jl?branch=master)

MDPTools is a package for solving discrete choice models often used in
economics. Currently, the package supports discrete and finite state spaces and
multinomial choice with linear in parameters utility functions. See [Aguirregabiria and Mira (2010)](http://www.sciencedirect.com/science/article/pii/S0304407609001985)
for a survey of the models this package is intended to handle. Standard methods
such as value function iterations (successive approximations), and policy/newton
iterations are implemented.

Note that for the extreme value shocks we use shocks with means normalized to zero.

![example](https://cloud.githubusercontent.com/assets/8431156/20754243/5360d216-b70a-11e6-906d-9eab9ed04d22.png)

The package is still in development. Things like continuous states, speed improvements
from providing additional structure (regenerative actions, ...), and relative value function
iterations are to be implemented in the future.

This package can be used on its own, but is meant to be a backend for estimation
packages yet to be announced: NestedFixedPoint.jl, NestedPseudoLikelihood.jl, MPEC.jl,
and a meta package to reexport this and the other packages called StructuralEstimation.jl.

This package is still in active development. Things will change.
