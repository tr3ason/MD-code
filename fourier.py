import matplotlib.pyplot as plt
import numpy as np
#Comment line
a = np.loadtxt("file.txt")

f = np.fft.rfft(a[:,1])
f=f/max(f)



fig, ax = plt.subplots(2, 1)
ax[0].plot(a[:,0],a[:,1])
ax[0].set_xlabel('Time')
ax[0].set_ylabel('Kinetic Energy')
n=len(f)
ax[1].plot(abs(f))
#ax[1].set_xlabel('Freq (Hz)')
#ax[1].set_ylabel('|Y(freq)|')
ax[1].set_ylim([0,1])
plt.show()
