#%%
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
import math

path='/storage/climatestor/PleioCEP/doensen/data/figs/'
mu = 0
variance = 1
sigma = math.sqrt(variance)
x = np.linspace(mu - 5*sigma, mu + 5*sigma, 100)
fig,ax=plt.subplots()
ax.plot(x, -stats.norm.pdf(x, mu, sigma),linewidth=4,color='k')
ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)
fig.savefig(path+'gaussian.png',dpi=300)

# %%
