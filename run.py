import subprocess
from multiprocessing import Pool


def run_process(kid):
    cmd = "scripts/kois-build-model -c butinah --submit {0}".format(kid)
    print("Running: {0}".format(cmd))
    subprocess.check_call(cmd, shell=True)


pool = Pool()
pool.map(run_process,
         [int(line.split()[0]) for line in open("targets/kois.txt")])
