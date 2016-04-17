#!/usr/bin/python

import subprocess

# Solve steady problem with L2 norm check
print("----------------------------------------------------")
print("| Steady Heat Equation - Gauss-Seidel with L2 Norm |")
print("----------------------------------------------------")
subprocess.call(["./main", "-p", "parametersL2.pot"])

# Solve steady problem with H1 norm check
print("----------------------------------------------------")
print("| Steady Heat Equation - Gauss-Seidel with H1 Norm |")
print("----------------------------------------------------")
subprocess.call(["./main", "-p", "parametersH1.pot"])

# Solve unsteady problem
print("---------------------------------------------")
print("| Unsteady Heat Equation - Thomas Algorithm |")
print("---------------------------------------------")
subprocess.call(["./main", "-p", "parametersUnsteady.pot"])
