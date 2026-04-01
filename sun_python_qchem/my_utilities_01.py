import time

def do_timing(fn, *args, **kwargs):
    t0 = time.time()
    fn(*args, **kwargs)
    t1 = time.time()
    t = t1 - t0
    print(f"{fn.__name__:25s} time = {t:10.5f} s")
