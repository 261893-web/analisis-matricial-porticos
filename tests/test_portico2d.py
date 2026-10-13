"""Validación contra soluciones cerradas de resistencia de materiales."""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from portico2d import Portico2D

E, A, I = 200e6, 0.01, 8e-5   # kN/m², m², m⁴


def test_voladizo_carga_puntual():
    # δ = PL³/3EI, θ = PL²/2EI
    P, L = 10.0, 3.0
    m = Portico2D()
    a, b = m.nudo(0, 0), m.nudo(L, 0)
    m.elemento(a, b, E, A, I)
    m.apoyo(a)
    m.carga(b, Fy=-P)
    m.resolver()
    ux, uy, rz = m.desplazamiento(b)
    assert uy == pytest.approx(-P * L**3 / (3 * E * I), rel=1e-9)
    assert rz == pytest.approx(-P * L**2 / (2 * E * I), rel=1e-9)
    assert m.reaccion(a)[2] == pytest.approx(P * L, rel=1e-9)


def test_viga_simple_carga_uniforme():
    # δmax = 5wL⁴/384EI, Mmax = wL²/8
    w, L = 12.0, 6.0
    m = Portico2D()
    a, c, b = m.nudo(0, 0), m.nudo(L / 2, 0), m.nudo(L, 0)
    m.elemento(a, c, E, A, I, w=-w)
    m.elemento(c, b, E, A, I, w=-w)
    m.apoyo(a, ux=True, uy=True, rz=False)
    m.apoyo(b, ux=False, uy=True, rz=False)
    m.resolver()
    assert m.desplazamiento(c)[1] == pytest.approx(-5 * w * L**4 / (384 * E * I), rel=1e-9)
    assert m.momento_en(0, L / 2) == pytest.approx(w * L**2 / 8, rel=1e-9)
    assert m.reaccion(a)[1] == pytest.approx(w * L / 2, rel=1e-9)


def test_viga_biempotrada():
    # Momento en el empotramiento = wL²/12
    w, L = 10.0, 5.0
    m = Portico2D()
    a, b = m.nudo(0, 0), m.nudo(L, 0)
    m.elemento(a, b, E, A, I, w=-w)
    m.apoyo(a)
    m.apoyo(b)
    m.resolver()
    assert abs(m.reaccion(a)[2]) == pytest.approx(w * L**2 / 12, rel=1e-9)
    assert m.momento_en(0, L / 2) == pytest.approx(w * L**2 / 24, rel=1e-9)


def test_equilibrio_global_portico():
    m = Portico2D()
    n0, n1, n2, n3 = m.nudo(0, 0), m.nudo(0, 4), m.nudo(6, 4), m.nudo(6, 0)
    m.elemento(n0, n1, E, A, I)
    m.elemento(n1, n2, E, A, I, w=-15)
    m.elemento(n2, n3, E, A, I)
    m.apoyo(n0)
    m.apoyo(n3)
    m.carga(n1, Fx=20)
    m.resolver()
    R = m.reaccion(n0) + m.reaccion(n3)
    assert R[0] == pytest.approx(-20, abs=1e-8)
    assert R[1] == pytest.approx(15 * 6, abs=1e-8)


def test_estructura_inestable():
    m = Portico2D()
    a, b = m.nudo(0, 0), m.nudo(4, 0)
    m.elemento(a, b, E, A, I)
    m.apoyo(a, ux=False, uy=True, rz=False)
    with pytest.raises(ValueError):
        m.resolver()
