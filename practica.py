#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sympy as sp
from sympy import symbols, Matrix
from sympy import init_printing
from sympy import simplify
from sympy import symbols, expand, factor

from pytc2.general import to_latex

from schemdraw import Drawing

# Ahora importamos las funciones de PyTC2

from pytc2.remociones import remover_polo_dc, remover_polo_jw , remover_polo_jw2
from pytc2.dibujar import display, dibujar_puerto_entrada, dibujar_funcion_exc_abajo,  dibujar_elemento_serie, dibujar_elemento_derivacion,  dibujar_tanque_derivacion, dibujar_tanque_RC_serie,  dibujar_espacio_derivacion, Capacitor, Resistor, ResistorIEC
from pytc2.dibujar import dibujar_Pi, dibujar_Tee, dibujar_lattice

from pytc2.general import print_latex, print_subtitle, a_equal_b_latex_s
from IPython.display import display,  Markdown

from pytc2.cuadripolos import Z2Tabcd_s, Y2Tabcd_s, Tabcd2Z_s, Tabcd2Y_s
from pytc2.cuadripolos import calc_MAI_impedance_ij, calc_MAI_vtransf_ij_mn, calc_MAI_ztransf_ij_mn
from pytc2.general import print_latex

# Activar la impresión en formato LaTeX
init_printing()
# Definir la variable simbólica s
s = sp.symbols('s',complex=True)

# Definir la función de transferencia H(s)
denominador = (s**2 + 2)*(s**2 + 5)
numerador = s *(3*s**2 + 7)
YY = numerador / denominador

print_subtitle('Admitancia de entrada a la red')

print_latex(a_equal_b_latex_s('Y(s)', YY))

##/////////////////////////////////////////////////////////////////////////////
##/////////////////////////////////////////////////////////////////////////////


# calculamos polos y ceros de Y(s)
num, den = YY.as_numer_denom()

roots_num = sp.solve(num, s, dict=True)
print(roots_num)

roots_den = sp.solve(den, s, dict=True)
print(roots_den)

# Restricción circuital: L2*C2 = 1 r/s
# remoción parcial en DC de 1/YY

Z2, Zc1 = remover_polo_dc(1/YY, omega_zero = 1 )

# Yc1 es la admitancia removida
# extraigo C1
C1 = 1/(s*Zc1)

print_latex(a_equal_b_latex_s('Z_1(s) = \\frac{k^p_0}{s}', Zc1))
print_latex(a_equal_b_latex_s('Z_2(s)', Z2))

##/////////////////////////////////////////////////////////////////////////////
##/////////////////////////////////////////////////////////////////////////////
Y4, Yt2, L2, C2 = remover_polo_jw2(1/Z2, isImpedance = False, omega = 1)

print_latex(a_equal_b_latex_s('Y_t3(s)', Yt2))
print_latex(a_equal_b_latex_s('Y_4(s)', Y4))

##/////////////////////////////////////////////////////////////////////////////
##/////////////////////////////////////////////////////////////////////////////

# quedaría solo un tanque en Y4, no especifico omega.
Y6, Yt3, L3, C3 = remover_polo_jw(Y4, isImpedance = False)

print_latex(a_equal_b_latex_s('Y_t5(s)', Yt3))
print_latex(a_equal_b_latex_s('Y_6(s)', Y6))

##/////////////////////////////////////////////////////////////////////////////
##/////////////////////////////////////////////////////////////////////////////
##/////////////////////////////////////////////////////////////////////////////
##/////////////////////////////////////////////////////////////////////////////

# Resolución simbólica
s = sp.symbols('s ', complex=True)

# Sea la siguiente función de excitación
YY = 3*s*(s**2+sp.Rational(7,3))/(s**2+2)/(s**2+5)

# Red ejemplo 1
d = dibujar_puerto_entrada(Drawing(unit=4),
                        voltage_lbl = ('+', '$V$', '-'), 
                        current_lbl = '$I$')

d = dibujar_funcion_exc_abajo(d, 
                 'Y',  
                 YY, 
                 hacia_salida = True,
                 k_gap_width = 0.5)

d = dibujar_elemento_serie(d, 'C', 'C1')

d = dibujar_tanque_derivacion(d, 'L2', 'C2')

d = dibujar_elemento_serie(d, 'L', 'L3')

d = dibujar_elemento_serie(d, 'C', 'C3')
display(d)

##/////////////////////////////////////////////////////////////////////////////
##/////////////////////////////////////////////////////////////////////////////

# Dibujamos la red resultante:

d = dibujar_puerto_entrada(Drawing(unit=4),
                        voltage_lbl = ('+', '$V$', '-'), 
                        current_lbl = '$I$')

d = dibujar_funcion_exc_abajo(d, 
                 'Y',  
                 YY, 
                 hacia_salida = True,
                 k_gap_width = 0.5)

d = dibujar_elemento_serie(d, 'C', C1)

d = dibujar_tanque_derivacion(d, L2, C2)

d = dibujar_elemento_serie(d, 'L', L3)

d = dibujar_elemento_serie(d, 'C', C3)

display(d)


##/////////////////////////////////////////////////////////////////////////////
##/////////////////////////////////////////////////////////////////////////////
##/////////////////////////////////////////////////////////////////////////////
##/////////////////////////////////////////////////////////////////////////////

print_subtitle('Calculos de verificacion')

Za, Zb, Zc = sp.symbols('Za, Zb, Zc', complex=True)

Ztee = sp.Matrix([[Za+Zb, Zb], [Zb, Zc + Zb]])

dibujar_Tee(Ztee)

Ztee_sym = sp.simplify(Ztee.subs(Za, 1/(s*C1)))
Ztee_sym = sp.simplify(Ztee_sym.subs(Zb, 1/(s*C2) + s*L2))
Ztee_sym = sp.simplify(Ztee_sym.subs(Zc, 1/(s*C3) + s*L3))

print_latex(a_equal_b_latex_s('Z_{tee}', Ztee_sym))

Ttee = Z2Tabcd_s(Ztee_sym)

Ytee = Tabcd2Y_s(Ttee)

print_latex(a_equal_b_latex_s('Y_{tee}', Ytee))

print_latex(a_equal_b_latex_s('Y_{11}', sp.factor(sp.simplify(sp.expand(Ytee[0])))) )

print_latex(a_equal_b_latex_s('Y_{21}', sp.factor(sp.simplify(sp.expand(Ytee[2])))) )

