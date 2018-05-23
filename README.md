# aboutDE


Let's gather here resources about Differential Equations to share among our
theoretical group.

## Presentations

- ODE *in silico*: https://github.com/KevCaz/ODEinsilico


## Code

### R

- packages required:

```r
install.packages(c("devtools", "deSolve", "odeint", "rmarkdown"))
devtools::install_github("inSileco/graphicsutils")
```

- scripts:

  - `R/RK4.R`: pure R (an naive) implementation of the [Runge–Kutta methods](https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods) (order 4) by @KevCaz.
  - `R/2Dsystems.R` Lotka-Volterra with the `deSolve` package by Torbjorg.
  - `R/vectorField.R`: a code to play with a vector field drawn for the linear system in 2 dimensions by @KevCaz.
  - `R/ODEwithR.Rmd`: let's build up a document to gather some thoughts!

### Julia

- [DifferentialEquations.jl](http://docs.juliadiffeq.org/latest/)

### Mathematica

- Mathematica: a side note to mention that you get it for free if you use a [raspberry pi](http://www.wolfram.com/raspberry-pi/)


- :video_camera: multimedia:

  - http://www.chaos-math.org/en: CHAOS A MATHEMATICAL ADVENTURE
