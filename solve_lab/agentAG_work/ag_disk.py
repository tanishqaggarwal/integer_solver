#!/usr/bin/env python3
"""AG: measure this box's actual random-read rate, to audit AB's
   'disk is not a way out: ~10^2/s seek rate vs ~10^8/s memory rate, a 2^20 slowdown'.

Method: write a 200 MB scratch file (inside AG's 300 MB budget, deleted at the end), open it
with O_DIRECT so the page cache is bypassed entirely, and time N random 4 KiB reads.
Also time N random reads from an in-RAM bytearray of the same size for the memory reference.
"""
import os, sys, time, random, ctypes, mmap

PATH = "/home/user/integer_solver/solve_lab/agentAG_work/_ag_seek.tmp"
SIZE = 200*1024*1024
BLK  = 4096
NREAD = 4000

def make():
    with open(PATH, "wb") as f:
        chunk = os.urandom(1 << 20)
        for _ in range(SIZE >> 20):
            f.write(chunk)
        f.flush(); os.fsync(f.fileno())

def bench_direct():
    O_DIRECT = getattr(os, "O_DIRECT", 0o40000)
    try:
        fd = os.open(PATH, os.O_RDONLY | O_DIRECT)
    except OSError as e:
        return None, "O_DIRECT unavailable: %s" % e
    # aligned buffer
    buf = mmap.mmap(-1, BLK)
    nblk = SIZE // BLK
    rnd = random.Random(20260807)
    offs = [rnd.randrange(nblk)*BLK for _ in range(NREAD)]
    t0 = time.perf_counter()
    for o in offs:
        os.preadv(fd, [buf], o)
    t1 = time.perf_counter()
    os.close(fd)
    return NREAD/(t1 - t0), None

def bench_ram():
    arr = bytearray(os.urandom(64*1024*1024))
    n = len(arr)//BLK
    rnd = random.Random(11)
    offs = [rnd.randrange(n)*BLK for _ in range(200000)]
    sink = 0
    t0 = time.perf_counter()
    for o in offs:
        sink ^= arr[o]
    t1 = time.perf_counter()
    return 200000/(t1 - t0), sink

if __name__ == "__main__":
    print("rotational flag reported by kernel for /dev/vda:",
          open("/sys/block/vda/queue/rotational").read().strip())
    make()
    try:
        iops, err = bench_direct()
        if err: print("  ", err)
        else:   print("MEASURED random 4KiB read rate, O_DIRECT (page cache bypassed): %.3e /s" % iops)
        rrate, _ = bench_ram()
        print("MEASURED random 4KiB touch rate from RAM (Python, so a LOWER bound):    %.3e /s" % rrate)
        if not err:
            import math
            print()
            print("AB claims disk ~1e2/s and memory ~1e8/s  ->  slowdown 2^%.1f" % math.log2(1e8/1e2))
            print("Measured ratio (RAM_python / disk) = 2^%.1f, and the true RAM rate is ~1e8-1e9/s,"
                  % math.log2(rrate/iops))
            print("so the honest slowdown against a C-speed hash table is 2^%.1f-2^%.1f."
                  % (math.log2(1e8/iops), math.log2(1e9/iops)))
    finally:
        if os.path.exists(PATH): os.remove(PATH)
        print("scratch file removed")
