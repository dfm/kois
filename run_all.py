import sqlite3
import subprocess


def run_process(kid):
    cmd = "scripts/kois-build-model {0} template.cfg --submit".format(kid)
    print("Running: {0}".format(cmd))
    subprocess.check_call(cmd, shell=True)


with sqlite3.connect("results/kois.db") as conn:
    c = conn.cursor()
    c.execute("select kepoi_name from kois where submitted is null")
    kids = map(lambda o: o[0], c.fetchall())

print(len(kids))
assert 0

map(run_process, kids)
