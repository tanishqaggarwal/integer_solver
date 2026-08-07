echo $$ > /home/user/integer_solver/solve_lab/agentT_work/job11.pid
export PYTHONDONTWRITEBYTECODE=1
cd /home/user/integer_solver/solve_lab/agentT_work
for s in 32 17 64; do
  echo "############ |S|=$s"
  /usr/bin/time -f "WALLCLOCK %e s" timeout 5400 python3 -u t_close2wj.py M$s $s 2>&1 | grep -v "^alignment\|^atoms with"
done
echo "ALLDONE"
