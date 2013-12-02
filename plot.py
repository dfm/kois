import os
import subprocess
# from multiprocessing import Pool


def run_process(kid):
    if os.path.exists("models/{0}/triangle.png".format(kid)):
        return None
    try:
        cmd = "scripts/kois-fetch-results {0}".format(kid)
        print("Running: {0}".format(cmd))
        subprocess.check_call(cmd, shell=True)

        cmd = "scripts/kois-plot-results models/{0} -b 20000".format(kid)
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
              [int(line.split()[0]) for line in open("targets/kois.txt")])
print("\n".join(map(str, [r for r in results if r is not None])))
