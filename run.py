import subprocess

def run_process(kid):
    cmd = "scripts/kois-build-model -c butinah --submit {0}".format(kid)
    print("Running: {0}".format(cmd))
    subprocess.check_call(cmd, shell=True)

map(run_process, [line.strip() for line in open("targets/kois.txt")])
