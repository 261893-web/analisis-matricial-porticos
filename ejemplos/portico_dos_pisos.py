"""
Ejemplo: pórtico de concreto armado de 2 pisos y 1 vano.
Cargas: gravedad repartida en vigas + carga lateral (sismo/viento simplificado).
Genera la deformada amplificada y el diagrama de momentos.
"""
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

RAIZ = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RAIZ / "src"))
from portico2d import Portico2D

# Concreto f'c = 210 kg/cm² -> E ≈ 15000*sqrt(210) kg/cm² ≈ 2.17e7 kN/m²
E = 2.17e7
col = dict(A=0.30 * 0.40, I=0.30 * 0.40**3 / 12)   # columna 30x40 cm
vig = dict(A=0.25 * 0.50, I=0.25 * 0.50**3 / 12)   # viga 25x50 cm

m = Portico2D()
L, h = 6.0, 3.2
n = [m.nudo(x, y) for y in (0, h, 2 * h) for x in (0, L)]
# n[0], n[1] base | n[2], n[3] piso 1 | n[4], n[5] piso 2

m.elemento(n[0], n[2], E, **col)
m.elemento(n[1], n[3], E, **col)
m.elemento(n[2], n[4], E, **col)
m.elemento(n[3], n[5], E, **col)
m.elemento(n[2], n[3], E, **vig, w=-25.0)   # kN/m
m.elemento(n[4], n[5], E, **vig, w=-18.0)

m.apoyo(n[0])
m.apoyo(n[1])
m.carga(n[2], Fx=15.0)
m.carga(n[4], Fx=25.0)

m.resolver()
print(m.resumen())

deriva = [(m.desplazamiento(n[2])[0]) / h,
          (m.desplazamiento(n[4])[0] - m.desplazamiento(n[2])[0]) / h]
for k, d in enumerate(deriva, 1):
    print(f"Deriva de entrepiso {k}: {d:.5f}  (elástica; para E.030 multiplicar por 0.75R antes de comparar con 0.007)")


# ---------------- gráficos ----------------
def deformada(m, idx, escala, npts=30):
    """Interpola la deformada con funciones de forma de Hermite."""
    e = m.elementos[idx]
    L, c, s = m._geometria(e)
    u = m._T(c, s) @ m.U[m._gdl(e)]
    xi = np.linspace(0, 1, npts)
    ux = (1 - xi) * u[0] + xi * u[3]
    N1, N2 = 1 - 3 * xi**2 + 2 * xi**3, L * (xi - 2 * xi**2 + xi**3)
    N3, N4 = 3 * xi**2 - 2 * xi**3, L * (-xi**2 + xi**3)
    uy = N1 * u[1] + N2 * u[2] + N3 * u[4] + N4 * u[5]
    x0, y0 = m.nudos[e["i"]]
    xl, yl = xi * L + escala * ux, escala * uy
    return x0 + c * xl - s * yl, y0 + s * xl + c * yl


fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 5.5))
for k, e in enumerate(m.elementos):
    (xi_, yi_), (xj_, yj_) = m.nudos[e["i"]], m.nudos[e["j"]]
    for ax in (ax1, ax2):
        ax.plot([xi_, xj_], [yi_, yj_], color="0.75", lw=2)
    ax1.plot(*deformada(m, k, escala=150), color="#1f5fa8", lw=2)

    # diagrama de momentos dibujado perpendicular al elemento
    L, c, s = m._geometria(e)
    x = np.linspace(0, L, 40)
    M = np.array([m.momento_en(k, xx) for xx in x])
    esc = 0.009
    px = xi_ + c * x + s * M * esc
    py = yi_ + s * x - c * M * esc
    ax2.fill(np.r_[xi_ + c * x, px[::-1]], np.r_[yi_ + s * x, py[::-1]],
             color="#c0392b", alpha=0.25)
    ax2.plot(px, py, color="#c0392b", lw=1.5)
    kmax = np.argmax(np.abs(M))
    ax2.annotate(f"{M[kmax]:.1f}", (px[kmax], py[kmax]), fontsize=8)

ax1.set_title("Deformada (amplificada x150)")
ax2.set_title("Momento flector [kN·m]")
for ax in (ax1, ax2):
    ax.set_aspect("equal")
    ax.axis("off")
plt.tight_layout()
plt.savefig(RAIZ / "figuras" / "portico_dos_pisos.png", dpi=150)
print("Figura guardada en figuras/portico_dos_pisos.png")
