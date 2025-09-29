# Black-Scholes-Pricer
Simple pricer using Black-Scholes, i will continue to update it, and add new features.
## What to expect from the program: 
-Greeks: Delta, Gamma, Vega, Theta (per year & per day), Rho  
-Handles special cases:
  - maturity now (T = 0) = intrinsic value
  - zero volatility (σ = 0) = parity lower bound
## How to use, user part
I added a new way of getting the parameters, using QtWidgets, this is really useful but it can creates some problems, si i had to add a lot of controls to make sure it will work in any cases
