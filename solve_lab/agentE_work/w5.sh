until [ $(wc -l < runs/scanA5.log) -ge 60 ]; do sleep 20; done
