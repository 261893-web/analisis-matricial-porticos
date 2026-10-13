"""
Análisis matricial de pórticos planos por el método de rigidez directa.

Cada nudo tiene 3 grados de libertad: (ux, uy, giro).
Cada elemento es una viga-columna de Euler-Bernoulli con deformación axial.

Unidades: el código no impone unidades, solo pide que sean consistentes
(en los ejemplos se usa kN y m).
"""

import numpy as np


class Portico2D:
    def __init__(self):
        self.nudos = []          # lista de (x, y)
        self.elementos = []      # lista de dicts con i, j, E, A, I, w
        self.apoyos = {}         # nudo -> (rx, ry, rz) con True = restringido
        self.cargas = {}         # nudo -> [Fx, Fy, M]

    # ---------- definición del modelo ----------
    def nudo(self, x, y):
        self.nudos.append((float(x), float(y)))
        return len(self.nudos) - 1

    def elemento(self, i, j, E, A, I, w=0.0):
        """w: carga distribuida uniforme perpendicular al eje (en ejes locales,
        positiva hacia +y local). Para una viga de izquierda a derecha,
        w negativa = carga de gravedad."""
        self.elementos.append(dict(i=i, j=j, E=E, A=A, I=I, w=w))
        return len(self.elementos) - 1

    def apoyo(self, n, ux=True, uy=True, rz=True):
        self.apoyos[n] = (ux, uy, rz)

    def carga(self, n, Fx=0.0, Fy=0.0, M=0.0):
        actual = self.cargas.get(n, [0.0, 0.0, 0.0])
        self.cargas[n] = [actual[0] + Fx, actual[1] + Fy, actual[2] + M]

    # ---------- matrices de elemento ----------
    def _geometria(self, e):
        xi, yi = self.nudos[e["i"]]
        xj, yj = self.nudos[e["j"]]
        L = np.hypot(xj - xi, yj - yi)
        c, s = (xj - xi) / L, (yj - yi) / L
        return L, c, s

    @staticmethod
    def _k_local(E, A, I, L):
        a = E * A / L
        b = 12 * E * I / L**3
        c = 6 * E * I / L**2
        d = 4 * E * I / L
        f = 2 * E * I / L
        return np.array([
            [ a,  0,  0, -a,  0,  0],
            [ 0,  b,  c,  0, -b,  c],
            [ 0,  c,  d,  0, -c,  f],
            [-a,  0,  0,  a,  0,  0],
            [ 0, -b, -c,  0,  b, -c],
            [ 0,  c,  f,  0, -c,  d],
        ])

    @staticmethod
    def _T(c, s):
        r = np.array([[c, s, 0], [-s, c, 0], [0, 0, 1]])
        T = np.zeros((6, 6))
        T[:3, :3] = r
        T[3:, 3:] = r
        return T

    @staticmethod
    def _empotramiento(w, L):
        """Fuerzas de empotramiento perfecto (locales) por carga uniforme w."""
        return np.array([0, w * L / 2, w * L**2 / 12, 0, w * L / 2, -w * L**2 / 12])

    def _gdl(self, e):
        i, j = e["i"], e["j"]
        return [3 * i, 3 * i + 1, 3 * i + 2, 3 * j, 3 * j + 1, 3 * j + 2]

    # ---------- solución ----------
    def resolver(self):
        n = 3 * len(self.nudos)
        K = np.zeros((n, n))
        F = np.zeros(n)

        for nd, (fx, fy, m) in self.cargas.items():
            F[3 * nd:3 * nd + 3] += [fx, fy, m]

        for e in self.elementos:
            L, c, s = self._geometria(e)
            T = self._T(c, s)
            kg = T.T @ self._k_local(e["E"], e["A"], e["I"], L) @ T
            g = self._gdl(e)
            K[np.ix_(g, g)] += kg
            # la carga distribuida entra como fuerzas nodales equivalentes
            F[g] += T.T @ self._empotramiento(e["w"], L)

        restringidos = [3 * nd + k for nd, r in self.apoyos.items()
                        for k in range(3) if r[k]]
        libres = np.setdiff1d(np.arange(n), restringidos)

        if np.linalg.matrix_rank(K[np.ix_(libres, libres)]) < len(libres):
            raise ValueError("Estructura inestable: revisa apoyos y conexiones.")

        U = np.zeros(n)
        U[libres] = np.linalg.solve(K[np.ix_(libres, libres)], F[libres])

        self.U = U
        self.reacciones = K @ U - F
        self.reacciones[libres] = 0.0
        self.fuerzas_internas = [self._fuerzas_elemento(e) for e in self.elementos]
        return U

    def _fuerzas_elemento(self, e):
        """Fuerzas en extremos del elemento en ejes locales:
        [N_i, V_i, M_i, N_j, V_j, M_j] (convención de la matriz de rigidez)."""
        L, c, s = self._geometria(e)
        T = self._T(c, s)
        u_local = T @ self.U[self._gdl(e)]
        return self._k_local(e["E"], e["A"], e["I"], L) @ u_local - self._empotramiento(e["w"], L)

    # ---------- utilidades de salida ----------
    def desplazamiento(self, nd):
        return self.U[3 * nd:3 * nd + 3]

    def reaccion(self, nd):
        return self.reacciones[3 * nd:3 * nd + 3]

    def momento_en(self, idx, x):
        """Momento flector a una distancia x del nudo i (convención de diseño:
        positivo = tracción en la fibra inferior local)."""
        e = self.elementos[idx]
        f = self.fuerzas_internas[idx]
        return -f[2] + f[1] * x + e["w"] * x**2 / 2

    def resumen(self):
        lineas = ["Nudo |      ux        uy       giro"]
        for k in range(len(self.nudos)):
            ux, uy, rz = self.desplazamiento(k)
            lineas.append(f"{k:4d} | {ux: .3e} {uy: .3e} {rz: .3e}")
        lineas.append("\nReacciones (Fx, Fy, M):")
        for k in self.apoyos:
            fx, fy, m = self.reaccion(k)
            lineas.append(f"  nudo {k}: {fx: 9.3f} {fy: 9.3f} {m: 9.3f}")
        return "\n".join(lineas)
