#!/usr/bin/python

import subprocess

# Run program with L2 norm check
subprocess.Popen(["./main", "-p", "parametersL2.pot"])

# Run program with H1 norm check
subprocess.Popen(["./main", "-p", "parametersH1.pot"])
