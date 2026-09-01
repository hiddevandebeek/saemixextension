# Publication code audit

Date: 31 August 2026

## Round 1: execution path

The supported estimator was traced from `gaussianCopulaFrem()` through the
saemix E- and M-step hooks. Twelve high-value regression tests established the
pre-refactor baseline. The public path reached 62 extension helpers but loaded
more than 200 extension functions.

## Round 2: mixed modules

`gaussianCopulaFrem.R` combined 773 lines of active density, conditioning and
sampling code with approximately 1,700 lines of collapse and margin-selection
experiments. The active block was retained. Correlation-coordinate helpers
were moved beside the other correlation helpers. The complete original file
was archived.

## Round 3: estimator boundary

The package runtime supported the publication score estimator and a historical
retained-particle common-Q estimator through the same global dispatch. The
runtime was replaced by a fixed-model score configuration. Common-Q pool,
optimizer, truncation, guard and delayed-activation branches were archived.
The saemix hooks now fail closed if a non-score runtime is encountered.

## Round 4: dependencies and interface

The C++ common-Q backend had no caller on the supported path. Its R wrappers,
C++ sources and build dependencies were removed from the package and archived.
Adaptive learning, direct natural-parameter margins and Theorem-8 prototypes
were also removed from package collation. Obsolete public documentation was
replaced by the fixed Gaussian-copula interface.

Repeated validation covered analytical scores, Fisher's identity, Gaussian
nesting, missing and categorical covariates, moving support, arbitrary
continuous eta margins, end-to-end fitting, likelihood evaluation and FFEM
translation. The clean staged package installs without native compilation.

## Remaining package-wide findings

The inherited `exportPattern()` exposes many original saemix helpers and all
alphabetically named extension internals, causing documentation warnings.
Several original saemix Rd signatures are also stale. These are package-wide
API maintenance issues rather than dependencies of the estimator. Removing
`exportPattern()` would change the established saemix API and was therefore not
done in this refactor. The separate undefined `msg` helper in the historical
discrete-VPC code was replaced by base `message()`.
