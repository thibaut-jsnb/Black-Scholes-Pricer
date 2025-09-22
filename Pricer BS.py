from math import log, sqrt, exp, erf, pi

S = 100.0     # Spot
K = 100.0     # Strike
T = 0.5       # Maturity (years)
r = 0.02      # risk-free rate
sigma = 0.25  # Volatility
q = 0.01      # Dividend (continuous)
opt_type = "call"  # "call" or "put"

# Mathematics tools
def N(x):
    return 0.5 * (1.0 + erf(x / sqrt(2.0)))

def n_pdf(x):
    return (1.0 / sqrt(2.0 * pi)) * exp(-0.5 * x * x)

# d1, d2 for Black–Scholes
def _d1_d2(S, K, T, r, sigma, q=0.0):
    if S <= 0 or K <= 0:
        raise ValueError("S and K have to be positive")
    if T <= 0 or sigma <= 0:
        return float('inf'), float('inf')
    vol_sqrtT = sigma * sqrt(T)
    d1 = (log(S / K) + (r - q + 0.5 * sigma * sigma) * T) / vol_sqrtT
    d2 = d1 - vol_sqrtT
    return d1, d2

# Black–Scholes price
def bs_price(S, K, T, r, sigma, q=0.0, option_type="call"):
    option_type = option_type.lower()
    if option_type not in ("call", "put"):
        raise ValueError("option_type must be named call or put.")
    if T <= 0:
        return max(S - K, 0.0) if option_type == "call" else max(K - S, 0.0)
    if sigma <= 0:
        if option_type == "call":
            return max(exp(-q * T) * S - exp(-r * T) * K, 0.0)
        else:
            return max(exp(-r * T) * K - exp(-q * T) * S, 0.0)
    d1, d2 = _d1_d2(S, K, T, r, sigma, q=q)
    if option_type == "call":     
        return exp(-q * T) * S * N(d1) - exp(-r * T) * K * N(d2)
    else:
        return exp(-r * T) * K * N(-d2) - exp(-q * T) * S * N(-d1)

# Greeks
def bs_greeks(S, K, T, r, sigma, q=0.0, option_type="call"):
    option_type = option_type.lower()
    if option_type not in ("call", "put"):
        raise ValueError("option_type must be named call or put.")
    if T <= 0 or sigma <= 0:
        price = bs_price(S, K, 0.0, r, max(sigma, 1e-12), q, option_type)
        return {
            "Price": price,
            "Delta": 0.0,
            "Gamma": 0.0,
            "Vega_per_vol": 0.0,
            "Vega_per_1pct": 0.0,
            "Theta_per_year": 0.0,
            "Theta_per_day": 0.0,
            "Rho_per_rate": 0.0,
            "Rho_per_bp": 0.0,
        }
    d1, d2 = _d1_d2(S, K, T, r, sigma, q=q)
    disc_r = exp(-r * T)
    disc_q = exp(-q * T)
    pdf_d1 = n_pdf(d1)

    # Delta
    if option_type == "call":
        delta = disc_q * N(d1)
    else:
        delta = disc_q * (N(d1) - 1.0)

    # Gamma
    gamma = (disc_q * pdf_d1) / (S * sigma * sqrt(T))

    # Vega
    vega_per_vol = disc_q * S * pdf_d1 * sqrt(T)
    vega_per_1pct = vega_per_vol / 100.0

    # Theta
    if option_type == "call":
        theta = (-disc_q * S * pdf_d1 * sigma / (2.0 * sqrt(T))
                 - r * disc_r * K * N(d2)
                 + q * disc_q * S * N(d1))
    else:
        theta = (-disc_q * S * pdf_d1 * sigma / (2.0 * sqrt(T))
                 + r * disc_r * K * N(-d2)
                 - q * disc_q * S * N(-d1))
    theta_per_year = theta
    theta_per_day = theta / 365.0  

    # Rho 
    if option_type == "call":
        rho = T * disc_r * K * N(d2)
    else:
        rho = -T * disc_r * K * N(-d2)
    rho_per_rate = rho
    rho_per_bp = rho / 10000.0  
    
    price = bs_price(S, K, T, r, sigma, q, option_type)
    return {
        "Price": price,
        "Delta": delta,
        "Gamma": gamma,
        "Vega_per_vol": vega_per_vol,
        "Vega_per_1pct": vega_per_1pct,
        "Theta_per_year": theta_per_year,
        "Theta_per_day": theta_per_day,
        "Rho_per_rate": rho_per_rate,
        "Rho_per_bp": rho_per_bp,
    }

print("Black Scholes European Option")
print(f"Parameters : S={S}, K={K}, T={T}, r={r}, sigma={sigma}, q={q}, type={opt_type}")

# Price + Greeks
price = bs_price(S, K, T, r, sigma, q, opt_type)
greeks = bs_greeks(S, K, T, r, sigma, q, opt_type)

print("\n Price and Greeks")
print(f"Price: {greeks['Price']:.6f}")
print(f"Delta: {greeks['Delta']:.6f}")
print(f"Gamma: {greeks['Gamma']:.6f}")
print(f"Vega (per 1.0 vol): {greeks['Vega_per_vol']:.6f}")
print(f"Vega (per +1% vol): {greeks['Vega_per_1pct']:.6f}")
print(f"Theta per year: {greeks['Theta_per_year']:.6f}")
print(f"Theta per day: {greeks['Theta_per_day']:.6f}")
print(f"Rho (per 1.0 rate): {greeks['Rho_per_rate']:.6f}")
print(f"Rho (per 1 bp): {greeks['Rho_per_bp']:.8f}")