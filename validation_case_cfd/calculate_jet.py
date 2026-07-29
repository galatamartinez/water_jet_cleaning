#!/usr/bin/env python3
"""
Jet Velocity and Area Calculator / Calculadora de Velocidad y Área del Chorro

Computes nozzle area, velocity, mass flow rate, and jet force for both:
1. Uniform Flow (Plug Flow): F_uniform = rho * A * U_avg^2 = m_dot * U_avg
2. Hagen-Poiseuille Flow (Laminar parabolic profile):
   - U_max = 2 * U_avg
   - Momentum flux coefficient beta = 4/3 ≈ 1.3333
   - F_poiseuille = (4/3) * rho * A * U_avg^2 = (4/3) * F_uniform
"""

import math
import sys
import argparse

def calculate_jet(flow_rate_lpm=60.0, diameter_mm=1.0, rho=1000.0):
    """
    Calculates jet hydraulic parameters and compares Uniform vs. Hagen-Poiseuille flows.
    
    Parameters:
    - flow_rate_lpm: Volumetric flow rate in L/min
    - diameter_mm: Nozzle inner diameter in mm
    - rho: Fluid density in kg/m³ (default: 1000 for water)
    
    Returns:
    - Dict with calculated parameters.
    """
    # Convert inputs to SI units
    Q_m3s = (flow_rate_lpm * 1e-3) / 60.0  # L/min to m³/s
    radius_m = (diameter_mm * 1e-3) / 2.0   # mm to m
    A_m2 = math.pi * (radius_m ** 2)         # Cross-sectional area in m²
    
    # Mean Velocity: U_avg = Q / A
    velocity_avg_ms = Q_m3s / A_m2
    
    # Mass flow rate: m_dot = rho * Q
    mass_flow_kgs = rho * Q_m3s
    
    # 1. Uniform Flow Momentum Force: F_uniform = m_dot * U_avg = rho * A * U_avg^2
    force_uniform_N = mass_flow_kgs * velocity_avg_ms
    
    # 2. Hagen-Poiseuille Flow (Parabolic Profile):
    # Velocity profile: u(r) = u_max * (1 - r^2 / R^2)
    # Peak center velocity: u_max = 2 * U_avg
    u_max_poiseuille_ms = 2.0 * velocity_avg_ms
    
    # Momentum flux: int(rho * u^2 dA) = (4/3) * rho * A * U_avg^2 = (4/3) * F_uniform
    beta_poiseuille = 4.0 / 3.0
    force_poiseuille_N = beta_poiseuille * force_uniform_N
    
    return {
        "flow_rate_lpm": flow_rate_lpm,
        "diameter_mm": diameter_mm,
        "rho_kg_m3": rho,
        "Q_m3s": Q_m3s,
        "A_m2": A_m2,
        "A_mm2": A_m2 * 1e6,
        "velocity_avg_ms": velocity_avg_ms,
        "mass_flow_kgs": mass_flow_kgs,
        "force_uniform_N": force_uniform_N,
        "u_max_poiseuille_ms": u_max_poiseuille_ms,
        "beta_poiseuille": beta_poiseuille,
        "force_poiseuille_N": force_poiseuille_N
    }

def print_summary(res):
    print("=" * 70)
    print(" CÁLCULO DE PARÁMETROS DEL CHORRO / JET PARAMETERS COMPARISON")
    print("=" * 70)
    print(f"  Caudal (Q):                   {res['flow_rate_lpm']:.2f} L/min ({res['Q_m3s']:.6e} m³/s)")
    print(f"  Diámetro de la tobera:        {res['diameter_mm']:.4f} mm")
    print(f"  Densidad del fluido (rho):    {res['rho_kg_m3']:.1f} kg/m³")
    print(f"  Área transversal (A):         {res['A_m2']:.6e} m² ({res['A_mm2']:.4f} mm²)")
    print(f"  Flujo másico (m_dot):         {res['mass_flow_kgs']:.4f} kg/s")
    print("-" * 70)
    print(" 1. FLUJO UNIFORME (Plug / Uniform Flow):")
    print(f"     Velocidad media (U_avg):   {res['velocity_avg_ms']:.2f} m/s")
    print(f"     Fuerza de impulso (F_unif): {res['force_uniform_N']:.2f} N  (100.0%)")
    print("-" * 70)
    print(" 2. FLUJO HAGEN-POISEUILLE (Parabólico / Laminar Flow):")
    print(f"     Velocidad máxima (u_max):   {res['u_max_poiseuille_ms']:.2f} m/s  (2.0 x U_avg)")
    print(f"     Factor de impulso (beta):   {res['beta_poiseuille']:.4f} (4/3)")
    print(f"     Fuerza de impulso (F_Pois): {res['force_poiseuille_N']:.2f} N  (133.3% vs Uniforme)")
    print("=" * 70)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculadora de velocidad y área del chorro")
    parser.add_argument("-q", "--flow-rate", type=float, default=4.8, help="Caudal en L/min (defecto: 60.0)")
    parser.add_argument("-d", "--diameter", type=float, default=1.0, help="Diámetro de tobera en mm (defecto: 1.0)")
    parser.add_argument("-r", "--rho", type=float, default=1000.0, help="Densidad del fluido en kg/m³ (defecto: 1000.0)")
    
    args = parser.parse_args()
    results = calculate_jet(flow_rate_lpm=args.flow_rate, diameter_mm=args.diameter, rho=args.rho)
    print_summary(results)
