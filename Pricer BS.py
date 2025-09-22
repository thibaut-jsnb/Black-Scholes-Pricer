# -*- coding: utf-8 -*-
"""
Pricer Black-Scholes (European Options) 

- Prix Call/Put (Black–Scholes, rendement de dividende continu facultatif q)
- Greeks : Delta, Gamma, Vega, Theta, Rho
- Vol implicite (méthode de bissection, robuste et sans dépendances)
- Mini démo de couverture delta (scénarios de variation du sous-jacent)
 
⚠️ Unit tests "visuels" à la fin dans le bloc __main__, à adapter à ton cas.
T = maturité en années (ex: 0.5 = 6 mois)
r = taux sans risque (continu) en décimal (ex: 0.02 = 2%)
q = rendement de dividende (continu) en décimal (ex: 0.01 = 1%)

A coller tel quel dans Spyder.
"""

from math import log, sqrt, exp, erf, pi

# =========================
# 1) Outils mathématiques
# =========================

def N(x):
    """
    Fonction de répartition de la loi normale standard : N(x)
    Implémentée via erf (fonction d'erreur). Évite SciPy.
    """
    return 0.5 * (1.0 + erf(x / sqrt(2.0)))

def n_pdf(x):
    """
    Densité de la loi normale standard : φ(x)
    """
    return (1.0 / sqrt(2.0 * pi)) * exp(-0.5 * x * x)

# =========================
# 2) d1, d2 de Black-Scholes
# =========================

def _d1_d2(S, K, T, r, sigma, q=0.0):
    """
    Calcule d1 et d2 classiques de Black–Scholes.
    - S : prix spot
    - K : strike
    - T : maturité en années
    - r : taux sans risque (continu)
    - sigma : volatilité (en décimal, ex: 0.2 pour 20%)
    - q : dividende (rendement continu), 0 par défaut

    Gestion de cas limites :
    - si T <= 0 : on est à maturité => pas de volatilité temporelle (traité ailleurs)
    - si sigma <= 0 : volatilité nulle => cas limite (traité ailleurs)
    """
    if S <= 0 or K <= 0:
        raise ValueError("S et K doivent être strictement positifs.")
    if T <= 0 or sigma <= 0:
        # On renvoie quelque chose (pas utilisé si on gère correctement les bornes)
        return float('inf'), float('inf')

    # standard formula
    # d1 = [ln(S/K) + (r - q + 0.5 sigma^2) T] / (sigma sqrt(T))
    # d2 = d1 - sigma sqrt(T)
    vol_sqrtT = sigma * sqrt(T)
    d1 = (log(S / K) + (r - q + 0.5 * sigma * sigma) * T) / vol_sqrtT
    d2 = d1 - vol_sqrtT
    return d1, d2

# Black Scholes Price

def bs_price(S, K, T, r, sigma, q=0.0, option_type="call"):
   
    option_type = option_type.lower()
    if option_type not in ("call", "put"):
        raise ValueError("option_type doit être 'call' ou 'put'.")

    # Instant Maturity Case
    if T <= 0:    
        if option_type == "call":
            return max(S - K, 0.0)
        else:
            return max(K - S, 0.0)

    # Sigma 0 case
    if sigma <= 0:
        if option_type == "call":
            return max(exp(-q * T) * S - exp(-r * T) * K, 0.0)
        else:
            return max(exp(-r * T) * K - exp(-q * T) * S, 0.0)

    d1, d2 = _d1_d2(S, K, T, r, sigma, q=q)

    if option_type == "call":
        # C = e^{-qT} S N(d1) - e^{-rT} K N(d2)
        return exp(-q * T) * S * N(d1) - exp(-r * T) * K * N(d2)
    else:
        # P = e^{-rT} K N(-d2) - e^{-qT} S N(-d1)
        return exp(-r * T) * K * N(-d2) - exp(-q * T) * S * N(-d1)

# =========================
# 4) Greeks analytiques
# =========================

def bs_greeks(S, K, T, r, sigma, q=0.0, option_type="call"):
    """
    Calcule les sensibilités (Greeks) analytiques.
    Retourne un dict : Delta, Gamma, Vega, Theta (par jour et par an), Rho.
    Vega est renvoyé par 'point de vol' (1.0 = 100% de vol) et par '1%' (divisé par 100).
    Theta renvoyé par an (convention continue) ET par jour (365).
    Rho par unité de taux (1.0 = +100% de taux), et par '1 bp' (divisé par 10000).
    """
    option_type = option_type.lower()
    if option_type not in ("call", "put"):
        raise ValueError("option_type doit être 'call' ou 'put'.")

    # Cas limites : même logique que pour le prix
    if T <= 0 or sigma <= 0:
        # Greeks mal définis au point T=0 ; on renvoie des zéros cohérents
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

    # --- Delta ---
    if option_type == "call":
        delta = disc_q * N(d1)
    else:
        delta = disc_q * (N(d1) - 1.0)

    # --- Gamma ---
    gamma = (disc_q * pdf_d1) / (S * sigma * sqrt(T))

    # --- Vega ---
    # Vega = ∂Price/∂sigma ; par "point de vol" (ex: 0.01 -> 1% vol = vega/100)
    vega_per_vol = disc_q * S * pdf_d1 * sqrt(T)
    vega_per_1pct = vega_per_vol / 100.0

    # --- Theta ---
    # Formules classiques (en convention "par an")
    if option_type == "call":
        theta = (-disc_q * S * pdf_d1 * sigma / (2.0 * sqrt(T))
                 - r * disc_r * K * N(d2)
                 + q * disc_q * S * N(d1))
    else:
        theta = (-disc_q * S * pdf_d1 * sigma / (2.0 * sqrt(T))
                 + r * disc_r * K * N(-d2)
                 - q * disc_q * S * N(-d1))

    theta_per_year = theta
    theta_per_day = theta / 365.0  # convention calendrier simple

    # --- Rho ---
    if option_type == "call":
        rho = T * disc_r * K * N(d2)
    else:
        rho = -T * disc_r * K * N(-d2)

    rho_per_rate = rho
    rho_per_bp = rho / 10000.0  # par "basis point" (1 bp = 0.01%)

    # Prix pour info
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

# =========================
# 5) Volatilité implicite
# =========================

def implied_vol_bisection(target_price, S, K, T, r, q=0.0, option_type="call",
                          vol_lower=1e-6, vol_upper=5.0, tol=1e-7, max_iter=200):
    """
    Trouve la volatilité implicite par bissection.
    - target_price : prix observé sur le marché
    - on cherche sigma tel que BS_price(sigma) ≈ target_price
    - vol_lower/vol_upper : bornes de recherche (par défaut [~0%, 500%])
    - tol : tolérance absolue sur l'erreur de prix
    - max_iter : itérations max

    Stratégie :
    1) Vérifie que f(lower) et f(upper) encadrent la racine (produit < 0).
    2) Bissection jusqu'à tol.
    3) Retourne sigma.
    """
    if T <= 0:
        raise ValueError("T doit être > 0 pour l'implied vol.")
    if target_price <= 0:
        raise ValueError("Le prix cible doit être positif.")

    def f(sig):
        return bs_price(S, K, T, r, sig, q, option_type) - target_price

    f_low = f(vol_lower)
    f_high = f(vol_upper)

    # Si pas d'encadrement, on tente d'élargir (optionnel)
    if f_low * f_high > 0:
        # Essai d'élargissement basique
        for factor in [2, 5, 10]:
            f_high = f(vol_upper * factor)
            if f_low * f_high <= 0:
                vol_upper *= factor
                break
        else:
            raise ValueError("Impossible d'encadrer la racine pour l'implied vol.")

    # Bissection
    a, b = vol_lower, vol_upper
    for _ in range(max_iter):
        m = 0.5 * (a + b)
        fm = f(m)
        if abs(fm) < tol:
            return m
        # On garde l'intervalle qui contient la racine
        if f_low * fm <= 0:
            b = m
            f_high = fm
        else:
            a = m
            f_low = fm
    # Si on sort ici, on renvoie la meilleure approx (milieu)
    return 0.5 * (a + b)

# =========================
# 6) Mini démo couverture Delta
# =========================

def delta_hedge_scenarios(S, K, T, r, sigma, q=0.0, option_type="call",
                          shifts_pct=(-0.1, -0.05, 0.0, 0.05, 0.1)):
    """
    Génère un petit tableau de scénarios de variation immédiate de S (ΔS),
    pour illustrer la couverture delta :
    - On prend la position "short option + Delta actions" au temps t.
    - On regarde la P&L instantanée si S saute à S*(1+shift), T constant (choc instantané).
    - NB : c'est pédagogique (convexité/gamma non couverte => résidu).

    Retourne une liste de dicts.
    """
    # Greeks initiaux
    g = bs_greeks(S, K, T, r, sigma, q, option_type)
    price0 = g["Price"]
    delta0 = g["Delta"]

    results = []
    for sh in shifts_pct:
        S1 = S * (1.0 + sh)                     # nouveau spot
        price1 = bs_price(S1, K, T, r, sigma, q, option_type)  # reprice option
        # Portefeuille delta-hedgé au temps t:
        #   P = -Option + Delta0 * S
        # Après le choc:
        #   P' = -Option(S1) + Delta0 * S1
        pnl = (-price1 + delta0 * S1) - (-price0 + delta0 * S)
        results.append({
            "Shift_%": sh * 100.0,
            "S_new": S1,
            "OptionPrice_new": price1,
            "Delta_initial": delta0,
            "Hedge_PnL": pnl
        })
    return results

# =========================
# 7) Démo / point d’entrée
# =========================

if __name__ == "__main__":
    # -------- Paramètres d'exemple (à modifier) --------
    S = 100.0     # Spot
    K = 100.0     # Strike
    T = 0.5       # Maturité (années) -> 0.5 = 6 mois
    r = 0.02      # Taux sans risque (continu) = 2%
    sigma = 0.25  # Volatilité = 25%
    q = 0.01      # Dividende continu = 1%
    opt_type = "call"  # "call" ou "put"

    print("=== Black–Scholes European Option ===")
    print(f"Paramètres: S={S}, K={K}, T={T}, r={r}, sigma={sigma}, q={q}, type={opt_type}")

    # Prix + Greeks
    price = bs_price(S, K, T, r, sigma, q, opt_type)
    greeks = bs_greeks(S, K, T, r, sigma, q, opt_type)

    print("\n-- Prix & Greeks --")
    print(f"Price: {greeks['Price']:.6f}")
    print(f"Delta: {greeks['Delta']:.6f}")
    print(f"Gamma: {greeks['Gamma']:.6f}")
    print(f"Vega (par 1.0 de vol): {greeks['Vega_per_vol']:.6f}")
    print(f"Vega (par +1% de vol): {greeks['Vega_per_1pct']:.6f}")
    print(f"Theta par an: {greeks['Theta_per_year']:.6f}")
    print(f"Theta par jour: {greeks['Theta_per_day']:.6f}")
    print(f"Rho (par 1.0 de taux): {greeks['Rho_per_rate']:.6f}")
    print(f"Rho (par 1 bp): {greeks['Rho_per_bp']:.8f}")

    # Vol implicite à partir d'un prix cible (ex: le prix théorique + 0.50)
    target = price + 0.50
    try:
        iv = implied_vol_bisection(target, S, K, T, r, q, opt_type)
        print("\n-- Volatilité implicite (à partir d'un prix cible) --")
        print(f"Target price: {target:.6f} -> Implied vol: {iv:.6f}")
    except ValueError as e:
        print(f"[ImpliedVol Error] {e}")

    # Démo couverture delta
    print("\n-- Scénarios de couverture delta (choc instantané sur S) --")
    table = delta_hedge_scenarios(S, K, T, r, sigma, q, opt_type,
                                  shifts_pct=(-0.1, -0.05, 0.0, 0.05, 0.1))
    # Affichage propre
    print(f"{'Shift(%)':>8} | {'S_new':>10} | {'Price_new':>12} | {'Delta0':>8} | {'Hedge P&L':>12}")
    print("-"*62)
    for row in table:
        print(f"{row['Shift_%']:8.2f} | {row['S_new']:10.4f} | {row['OptionPrice_new']:12.4f} | "
              f"{row['Delta_initial']:8.4f} | {row['Hedge_PnL']:12.4f}")

    # Note pédagogique :
    # - La P&L de la couverture delta n'est pas nulle lorsque le choc est grand,
    #   car tu n'as pas couvert la Γ (gamma). Pour de "petits" dS, la P&L se rapproche de 0.


