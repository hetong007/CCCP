CCCP
=======

This is a repository implementing CCCP algortihm.

Run `data.R` on simulated data.

`lasso.cd`: Lasso with Coordinate Descent
`lasso.cd.cccp`: Algorithm from *Calibrating nonconvex penalized regression in ultra-high dimension*
`lasso.cd.ncv`: Algorithm from *Coordinate descent algorithms for nonconvex penalized regression, with applications to biological feature selection*

`HBIC`: For a CCCP estimation, calculate `HBIC`, `MSE`, `SCAD`, `TP`, `FP` on a sequence of `lambda`.

## Results from 100 simulations

Average MSE/TP/FP, lambda selected by HBIC/RMSE:

|      | HBIC/RMSE| MSE.bhat|   TP|    FP|
|:-----|---------:|--------:|----:|-----:|
|CCCP  |  1.974429| 2.641697| 2.65|  0.29|
|NCV   |  2.130215| 2.033658| 2.86|  0.43|
|Lasso |  2.229269| 1.084010| 3.00| 31.98|

Average MSE/TP/FP, lambda selected by MSE (given true beta):

|      | HBIC/RMSE|  MSE.bhat|   TP|    FP|
|:-----|---------:|---------:|----:|-----:|
|CCCP  |  4.937690| 0.5829823| 2.99| 32.35|
|NCV   |  4.982644| 0.6404177| 3.00| 31.83|
|Lasso |  2.294679| 0.8897132| 3.00| 31.98|
