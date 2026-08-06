#!/bin/sh
S=/home/user/integer_solver/solve_lab/s10
python3 $S/bl_scan.py canon c_pair_cone pair - cone      > $S/bl_scan_c_pair_cone.log 2>&1
python3 $S/bl_scan.py canon c_b00 single 2081 all        > $S/bl_scan_c_b00.log 2>&1
python3 $S/bl_scan.py canon c_b01 single 2081,4287 all   > $S/bl_scan_c_b01.log 2>&1
