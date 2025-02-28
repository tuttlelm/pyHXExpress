
import threading
import subprocess

def run_script(script_name):
    subprocess.run(["python3", script_name])

scripts = {'D109H': '/data/tuttle/HDX-MS/Chris_B5_disease/SpecExport/batchrun_2pop_15Aug2024/D109H_run/run_D109H_15Aug2024.py', 'R120G': '/data/tuttle/HDX-MS/Chris_B5_disease/SpecExport/batchrun_2pop_15Aug2024/R120G_run/run_R120G_15Aug2024.py', 'WT': '/data/tuttle/HDX-MS/Chris_B5_disease/SpecExport/batchrun_2pop_15Aug2024/WT_run/run_WT_15Aug2024.py'}

if __name__ == "__main__":
    script_thread = {}
    for k,v in scripts.items():
        script_thread[k] = threading.Thread(target=run_script, args=[v])
        #print("value",v)

    for k,v in scripts.items():
        script_thread[k].start()

    for k,v in scripts.items():
        script_thread[k].join()

    print("Batch scripts have finished executing.")
