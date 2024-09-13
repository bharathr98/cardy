import numpy as np
from matplotlib import pyplot as plt
import sys

exec(open(sys.argv[1], "rb").read())

data = np.array(data)
print(f'Number of samples - {len(data.flatten())}')

plt.hist(data.flatten(),bins=200);
plt.savefig(sys.argv[2])