import sqlite3
import subprocess


def run_process(kid):
    print("STARTING {0}".format(kid))
    try:
        cmd = ("scripts/kois-plot-results {0} {1} --triangle"
               .format(kid, "template.cfg"))
        print("Running: {0}".format(cmd))
        subprocess.check_call(cmd, shell=True)

    except Exception as e:
        print("{0} FAILED".format(kid))
        print(e)

    else:
        print("finished {0}".format(kid))
        return kid


with sqlite3.connect("results/kois.db") as conn:
    c = conn.cursor()
    c.execute("select kepoi_name from kois where submitted is not null")
    kids = map(lambda o: o[0], c.fetchall())

results = map(run_process, kids)
