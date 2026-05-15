#!/bin/bash
# Invoke driver.py with the Python interpreter from the dakota-env conda
# environment, which has both hydroravens and dakota.interfacing installed.
/home/awickert/anaconda3/envs/dakota-env/bin/python driver.py "$@"
