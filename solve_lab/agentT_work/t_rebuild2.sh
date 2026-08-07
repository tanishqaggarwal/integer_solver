set -e
set -o pipefail   # T37: a pipeline reports its LAST command, so `python3 X.py | tail -3`
                   # made `set -e` a no-op -- a crashing stage still reached "REBUILD DONE".
export PYTHONDONTWRITEBYTECODE=1
ML=/home/user/integer_solver/solve_lab/agentT_work/mirror/L
cd $ML
echo "=== buildall"; python3 -u buildall.py 2>&1 | tail -4
echo "=== calib2";   python3 -u calib2.py   2>&1 | tail -5
echo "=== slopes";   python3 -u slopes.py   2>&1 | tail -3
echo "=== DONE"; ls -la $ML/*.pkl
