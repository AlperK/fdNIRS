import time
from pathlib import Path
import numpy as np
from datetime import datetime


loc = Path('2023-11-16', 'DUAL-SLOPE-830-3', '1', 'times.txt')
timeStamps = np.loadtxt(fname=loc)

deltas = np.diff(timeStamps/8)
print(deltas)
# tic = datetime.timestamp(datetime.now())
# time.sleep(1)
# toc = datetime.timestamp(datetime.now())
# print(toc - tic)
