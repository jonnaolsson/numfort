# Source code

The original Fortran 77 source code can be downloaded from
[https://www.zhangzk.net/software.html](https://www.zhangzk.net/software.html).

## Numfort interface

- Optional `INFO` reports `NEWUOA_CONVERGED`, `NEWUOA_MAXFUN`,
  `NEWUOA_NUMERICAL_ERROR`, or `NEWUOA_INVALID_INPUT`.
- `NF` counts objective calls, including evaluations that return nonfinite values.
- Optional `FX` returns the residual already evaluated at the returned `X`.
  No extra objective call is needed. If the first evaluation fails, `FX` contains
  that residual; if input is invalid, `FX` is undefined.
- The best point and residual are retained during partial initialization and
  restored together on budget exhaustion or numerical failure.
- A trust-region step with a positive predicted change and length at least
  half the radius returns `NEWUOA_NUMERICAL_ERROR`. It does not terminate the
  caller or evaluate the rejected step. `minimize_dfls` reports
  `NF_STATUS_INVALID_STATE` with `success=false`.
