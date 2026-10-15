"""Wall time and peak memory of one run vs grid side, settled set, 16 C, single process."""
import os, resource, sys, time

from sweep_minimal import run_one

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.abspath(os.path.join(HERE, "..", "input_minimal_e1500"))
side, days = int(sys.argv[1]), int(sys.argv[2])
t0 = time.perf_counter()
run_one(SRC, {}, side, days, 0, temp=16)
dt = time.perf_counter() - t0
rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6  # macOS: bytes -> MB
print(f"side {side:3d} ({side*side:5d} cells) days {days}: {dt:7.1f} s total, {dt/days*1000:8.1f} ms/day, peak RSS {rss:6.0f} MB")
