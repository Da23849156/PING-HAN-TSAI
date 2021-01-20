import itertools
import numpy as np
import matplotlib as mpl
from sklearn.mixture import BayesianGaussianMixture
from collections import Counter
from scipy import linalg
import matplotlib.pyplot as plt

import contextlib   # context manager for printoptions later
@contextlib.contextmanager
def printoptions(*args, **kwargs):
    original = np.get_printoptions()
    np.set_printoptions(*args, **kwargs)
    try:
        yield
    finally: 
        np.set_printoptions(**original)

# scikit-learn.org/stable/auto_examples/mixture/plot_gmm.thml

color_iter = itertools.cycle(['navy', 'c', 'cornflowerblue', 'gold',
                              'darkorange'])
def plot_results(X, Y_, means, covariances, index, title):
    splot = plt.subplot(2, 1, 1 + index)
    for i, (mean, covar, color) in enumerate(zip(
            means, covariances, color_iter)):
        v, w = linalg.eigh(covar)
        v = 2. * np.sqrt(2.) * np.sqrt(v)
        u = w[0] / linalg.norm(w[0])
        # as the DP will not use every component it has access to
        # unless it needs it, we shouldn't plot the redundant
        # components.
        if not np.any(Y_ == i):
            continue
        plt.scatter(X[Y_ == i, 0], X[Y_ == i, 1], .8, color=color)

        # Plot an ellipse to show the Gaussian component
        angle = np.arctan(u[1] / u[0])
        angle = 180. * angle / np.pi  # convert to degrees
        ell = mpl.patches.Ellipse(mean, v[0], v[1], 180. + angle, color=color)
        ell.set_clip_box(splot.bbox)
        ell.set_alpha(0.5)
        splot.add_artist(ell)

    plt.xlim(-9., 5.)
    plt.ylim(-3., 6.)
    plt.xticks(())
    plt.yticks(())
    plt.title(title)


Y = np.genfromtxt('BGM_hypo4_data.csv', delimiter=',')

# rng = np.random.RandomState(random_state)
Y.shape

random_seed = 27132
n_components = 10
rng = np.random.RandomState(seed = random_seed)   # fix a seed

DProc = BayesianGaussianMixture(n_components=n_components,
  weight_concentration_prior_type="dirichlet_process",
  weight_concentration_prior=1e-1,
  n_init = 7,
  init_params='kmeans',      # default 'kmeans'
  random_state=random_seed   # if int, then taken as random seed
  ).fit(Y)   # random_state=random_state

results = DProc.predict(Y)
probs = DProc.predict_proba(Y)
res_prob = np.column_stack((probs, results))
# res_prob = np.around(res_prob, decimals = 3)  # around() for arrays

# context manager controls precision within the next block of print commands
with printoptions(precision=3, suppress=True):
    print(results)
    print("\nposterior prob:\n", probs)
    print("\nmean:\n", DProc.means_)
    print("\ncovariances\n", DProc.covariances_)
    print("\nweights", DProc.weights_)
    print("\nCount the clusters\n")
    print( Counter(results).keys() )     # equals to list(set(words))
    print( Counter(results).values() )   # count freq of the elements

np.savetxt('BGM_hypo4_Out.csv', results, fmt = '%1.1f', delimiter=',')
plot_results(Y, DProc.predict(Y), DProc.means_, DProc.covariances_, 1,
             'Bayesian Gaussian Mixture with a Dirichlet process prior')
plt.show()
