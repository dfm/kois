import os
import subprocess
# from multiprocessing import Pool


def run_process(kid):
    try:
        cmd = ("scripts/kois-plot-results {0} --burnin 20000 "
               "--triangle").format(kid)
        print("Running: {0}".format(cmd))
        subprocess.check_call(cmd, shell=True)

    except Exception as e:
        print("{0} FAILED".format(kid))
        print(e)

    else:
        print("FINISHED {0}".format(kid))
        return kid


# pool = Pool()
results = map(run_process,
              [float(line.strip()) for line in open("targets/kois.txt")])
print("\n".join(map(str, [r for r in results if r is not None])))
