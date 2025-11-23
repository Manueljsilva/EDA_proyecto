import math

def normal_cdf(x):
    # CDF de N(0,1) usando erf
    return 0.5 * (1.0 + math.erf(x / math.sqrt(2.0)))

def prob_colision(tau, w):
    # p(tau; w) = ∫_{-w/(2*tau)}^{w/(2*tau)} f(t) dt = CDF(b) - CDF(a)
    a = -w / (2.0 * tau)
    b =  w / (2.0 * tau)
    return normal_cdf(b) - normal_cdf(a)

def calcular_K_L(n, C, t, redondeo="ceil"):
    w0 = 4 * (C ** 2)

    p1 = prob_colision(1.0, w0)      # puntos cercanos
    p2 = prob_colision(C, w0)        # puntos lejanos
    rho = math.log(1.0 / p1) / math.log(1.0 / p2)

    base = 1.0 / p2
    objetivo = n / t                # según el paper (evita el error de t/n)
    K_real = math.log(objetivo) / math.log(base)

    if redondeo == "ceil":
        K = math.ceil(K_real)
    elif redondeo == "round":
        K = round(K_real)
    else:
        K = math.floor(K_real)

    L_real = (n / t) ** rho
    L = math.ceil(L_real)

    return {
        "w0": w0,
        "p1": p1,
        "p2": p2,
        "rho": rho,
        "K_real": K_real,
        "K": K,
        "L_real": L_real,
        "L": L,
    }

if __name__ == "__main__":
    n = 10000   # número de puntos a indexar
    C = 1.2   # ratio de aproximación
    t = 3000   # parámetro t

    res = calcular_K_L(n, C, t)
    for k, v in res.items():
        print(f"{k}: {v}")
