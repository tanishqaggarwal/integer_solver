#!/bin/sh
S=/home/user/integer_solver/solve_lab/s10
python3 $S/bl_scan.py f2 f2_b11 single 4287 all          > $S/bl_scan_f2_b11.log 2>&1
python3 $S/bl_scan.py f2 f2_b00 single 2081 all          > $S/bl_scan_f2_b00.log 2>&1
python3 $S/bl_scan.py f2 f2_b01 single 2081,4287 all     > $S/bl_scan_f2_b01.log 2>&1
python3 $S/bl_scan.py f2 f2_pair_pins pair - conepins    > $S/bl_scan_f2_pair_pins.log 2>&1
