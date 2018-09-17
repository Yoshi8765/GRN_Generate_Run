import matplotlib.pyplot as plt
import numpy as np

std = 0.05
mean = 0

errors = [np.random.normal(loc=mean, scale=std)]
errors_std = [errors[0]]
for i in range(100):
    next_error = np.random.normal(loc=mean, scale=std)
    errors.append(errors[i] + next_error)
    errors_std.append(next_error)

plt.subplot(311)
plt.plot(errors)

plt.subplot(312)
plt.plot(errors_std)
plt.show()
    
plt.subplot(313)

for i in range(1,len(errors)-1):
