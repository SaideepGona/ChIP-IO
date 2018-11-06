import os
import subprocess

DETACHED_PROCESS = 0x00000008

for x in range(10):
    print(x)

# os.system("python3 build_db.py")
# pid = subprocess.Popen("python3 build_db.py", shell=True, creationflags=DETACHED_PROCESS).pid
# os.spawnl(os.P_NOWAIT, "python3", "python3", "build_db.py")

for x in range(10):
    print(x)