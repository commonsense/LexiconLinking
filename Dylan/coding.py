import scipy
import numpy as np


x = np.eye(3)


print len(np.transpose(np.nonzero(x)))
print len(x == 0)
