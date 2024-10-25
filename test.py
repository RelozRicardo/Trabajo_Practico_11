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

from pytc2.remociones import remover_polo_dc, remover_polo_jw , remover_polo_jw2 ,remover_polo_infinito , remover_polo_dc2 , remover_polo_infinito2
from pytc2.dibujar import display, dibujar_tanque_serie, dibujar_puerto_entrada, dibujar_funcion_exc_abajo,  dibujar_elemento_serie, dibujar_elemento_derivacion,  dibujar_tanque_derivacion, dibujar_tanque_RC_serie,  dibujar_espacio_derivacion, Capacitor, Resistor, ResistorIEC
from pytc2.dibujar import dibujar_Pi, dibujar_Tee, dibujar_lattice

from pytc2.general import print_latex, print_subtitle, a_equal_b_latex_s
from IPython.display import display,  Markdown

from pytc2.cuadripolos import Z2Tabcd_s, Y2Tabcd_s, Tabcd2Z_s, Tabcd2Y_s
from pytc2.cuadripolos import calc_MAI_impedance_ij, calc_MAI_vtransf_ij_mn, calc_MAI_ztransf_ij_mn
from pytc2.general import print_latex

import scipy.signal as sig
from pytc2.sistemas_lineales import analyze_sys

# Activar la impresión en formato LaTeX
init_printing()
# Definir la variable simbólica s
s = sp.symbols('s',complex=True)

# Definir la función de transferencia H(s)
numerador = (s**2 + 1)*(s**2 + 3)
denominador = s *(s**2 + 2.5)
Z11 = numerador / denominador

print('Impedancia de entrada a la red')
print_latex(a_equal_b_latex_s('Z(s)', Z11))

##/////////////////////////////////////////////////////////////////////////////
##/////////////////////////////////////////////////////////////////////////////

Z2, Zt1 = remover_polo_infinito(Z11, omega_zero=3)
L1 = Zt1/s
L1 = sp.nsimplify(L1)
Z2 = sp.nsimplify(sp.factor(Z2))
Zt1 = sp.nsimplify(sp.factor(Zt1))

print_latex(a_equal_b_latex_s('Zt1(s)', Zt1))
print_latex(a_equal_b_latex_s('Z2(s)', Z2))

##/////////////////////////////////////////////////////////////////////////////
##/////////////////////////////////////////////////////////////////////////////

Y4, Yt3 , L2, C2 = remover_polo_jw(1/Z2, omega=3, isImpedance=False)
L2 = sp.nsimplify(L2)
C2 = sp.nsimplify(C2)
Y4 = sp.nsimplify(sp.factor(Y4))
Yt3 = sp.nsimplify(sp.factor(Yt3))

print_latex(a_equal_b_latex_s('Yt3(s)', Yt3))
print_latex(a_equal_b_latex_s('Y4(s)', Y4))

##/////////////////////////////////////////////////////////////////////////////
##/////////////////////////////////////////////////////////////////////////////

Z6, Zt5 = remover_polo_infinito2(1/Y4, omega_zero=2)
L3 = Zt5/s
L3 = sp.nsimplify(L3)
Zt5 = sp.nsimplify(sp.factor(Zt5))

print_latex(a_equal_b_latex_s('Zt5(s)', Zt5))
print_latex(a_equal_b_latex_s('Z6(s)', Z6))

##/////////////////////////////////////////////////////////////////////////////
##/////////////////////////////////////////////////////////////////////////////

Y8, Yt7 , L4, C4 = remover_polo_jw2(1/Z6, omega=2, isImpedance=False)
L4 = sp.nsimplify(L4)
C4 = sp.nsimplify(C4)
Y8 = sp.nsimplify(sp.factor(Y8))
Yt7 = sp.nsimplify(sp.factor(Yt7))

print_latex(a_equal_b_latex_s('Yt7(s)', Yt7))
print_latex(a_equal_b_latex_s('Y8(s)', Y8))


##/////////////////////////////////////////////////////////////////////////////
##/////////////////////////////////////////////////////////////////////////////
##/////////////////////////////////////////////////////////////////////////////
##/////////////////////////////////////////////////////////////////////////////

# Dibujamos la red resultante:

d = dibujar_puerto_entrada(Drawing(unit=4),
                        voltage_lbl = ('+', '$V$', '-'), 
                        current_lbl = '$I$')

d = dibujar_funcion_exc_abajo(d, 
                 'Z',  
                 Z11, 
                 hacia_salida = True,
                 k_gap_width = 0.5)

d = dibujar_elemento_serie(d, 'L', L1)

d = dibujar_tanque_derivacion(d, L2, C2)

d = dibujar_elemento_serie(d, 'L', L3)

d = dibujar_tanque_derivacion(d,  L4, C4)

display(d)

##/////////////////////////////////////////////////////////////////////////////
##/////////////////////////////////////////////////////////////////////////////
##/////////////////////////////////////////////////////////////////////////////
##/////////////////////////////////////////////////////////////////////////////

print('Calculos de verificacion')

Za1,Zb1 = sp.symbols('Za1 , Zb1', complex=True)
Za2,Zb2 = sp.symbols('Za2 , Zb2', complex=True)


T1 = sp.Matrix([[ 1 + Za1/Zb1 , Za1 ], [1/Zb1 , 1]])

T2 = sp.Matrix([[ 1 + Za2/Zb2 , Za2 ], [1/Zb2 , 1]])

print_latex(a_equal_b_latex_s('T_{1}', T1))
print_latex(a_equal_b_latex_s('T_{2}', T2))

TT = T1 * T2;

print_latex(a_equal_b_latex_s('T_{T}', TT))

TT_sym = sp.simplify(TT.subs(Za1, s*L1) ) 
TT_sym = sp.simplify(TT_sym.subs(Zb1, 1/(s*C2) + L2*s ))
TT_sym = sp.simplify(TT_sym.subs(Za2, s*L3 ) )
TT_sym = sp.simplify(TT_sym.subs(Zb2, 1/(s*C4) + L4*s ))

print_latex(a_equal_b_latex_s('T_{T}', TT_sym))


Tt_ori = 1/ TT_sym[0];
num , den = Tt_ori.as_numer_denom()


print_latex(a_equal_b_latex_s('TtEJ_{11}', sp.factor(sp.simplify(sp.expand(Tt_ori)))) )

##/////////////////////////////////////////////////////////////////////////////
##/////////////////////////////////////////////////////////////////////////////
##/////////////////////////////////////////////////////////////////////////////
##/////////////////////////////////////////////////////////////////////////////

V2V1 = sig.TransferFunction( [1, 0, 11, 0, 18], [1 , 0, 4, 0, 3] )

analyze_sys(V2V1)

