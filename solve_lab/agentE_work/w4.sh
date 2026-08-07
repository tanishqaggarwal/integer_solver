until [ $(wc -l < runs/scanA4.log) -ge 30 ]; do sleep 20; done
