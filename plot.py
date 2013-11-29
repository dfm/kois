import os
import subprocess
from multiprocessing import Pool


def run_process(kid):
    if os.path.exists("models/{0}/mcmc.txt".format(kid)):
        return None
    cmd = "scripts/kois-fetch-results {0}".format(kid)
    print("Running: {0}".format(cmd))
    subprocess.check_call(cmd, shell=True)

    cmd = "scripts/kois-plot-results {0} -b 20000".format(kid)
    print("Running: {0}".format(cmd))
    subprocess.check_call(cmd, shell=True)

    print("FINISHED {0}".format(kid))
    return kid


pool = Pool()
results = pool.map(run_process,
                   [int(line.split()[0]) for line in open("targets/kois.txt")])
print("\n".join(map(str, [r for r in results if r is not None])))
