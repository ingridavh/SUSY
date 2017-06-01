import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

mg = np.array([800, 1200, 1600, 2000])
mq = np.array([1000, 1500, 2000])
s_sqrt = np.linspace(3000, 13000, 11)

#mg=800, mq=1000 start at s=3000, 4000, 5000, 6000, 7000, 8000, 9000, 10 000, 11 000, 12 000, 13 000
LO = np.array([0.274E-05, 0.347E-03, 0.339E-02, 0.133E-01, 0.332E-01, 0.647E-01, 0.108, 0.162, 0.226, 0.299, 0.379])
NLO = np.array([0.368E-05, 0.428E-03, 0.414E-02, 0.162E-01, 0.408E-01, 0.798E-01, 0.134, 0.201, 0.282, 0.374, 0.476])







#mg=2000, mq=800 start at 4000, 5000, 6000, 7000, 8000, 9000, 10 000, 11 000
#12 000, 13 000
s_sqrt2 = np.linspace(4000, 13000, 10)
LO2 = np.array([0.250E-02, 0.124E-01, 0.336E-01, 0.672E-01, 0.113, 0.169, 0.234, 0.307, 0.386, 0.472])
NLO2 = np.array([0.243E-02, 0.120E-01, 0.328E-01, 0.662E-01, 0.112, 0.168, 0.235, 0.311, 0.394, 0.485])

sns.set_context("notebook", font_scale=1.5, rc={"lines.linewidth": 2.5})

plt.plot(s_sqrt2, LO2, s_sqrt2, NLO2)
plt.title('Cross section Prospino', size='large')
plt.legend(['LO', 'NLO'], loc='upper left', fontsize='large')
plt.xlabel('sqrt(s) [GeV]', size='large')
plt.ylabel('sigma', size='large')
plt.savefig('prospino_2000_800_tester.pdf')
plt.show()

plt.plot(s_sqrt2, NLO2-LO2)
plt.title('Cross section Prospino', size='large')
plt.legend(['NLO-LO'], loc='upper left', fontsize='large')
plt.xlabel('sqrt(s) [GeV]', size='large')
plt.ylabel('sigma', size='large')
plt.savefig('prospino_2000_800_tester_diff.pdf')
plt.show()

#plt.plot(s_sqrt, NLO-LO)
#plt.title('Difference in cross section, LO and NLO')
#plt.xlabel('sqrt(s) [GeV]')
#plt.ylabel('d sigma (NLO-LO)')
#plt.show()
