for a7 in 24601 112 24734 24774 33095; do
  for a34 in 2081 22562 28713 20564; do
    echo -n "a7=$a7 a34=$a34: "
    timeout 130 python3 fc_param.py $a7 $a34 2>&1 | grep FINAL
  done
done
