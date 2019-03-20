import numpy as np

from .r_function import RFunction

_SplatSimulate = RFunction(
    name="splat", setup="""
        library(splatter)
        library(scater)
        library(magrittr)
    """,
    args="""
        method='paths',
        nBatches=1, batchCells=100,
        nGenes=10000,
        batch_facLoc=0.1, batch_facScale=0.1,
        mean_rate=0.3, mean_shape=0.6,
        lib_loc=11, lib_scale=0.2, lib_norm=False,
        out_prob=0.05,
        out_facLoc=4, out_facScale=0.5,
        de_prob=0.1, de_downProb=0.1,
        de_facLoc=0.1, de_facScale=0.4,
        bcv_common=0.1, bcv_df=60,
        dropout_type='none', dropout_prob=0.5,
        dropout_mid=0, dropout_shape=-1,
        group_prob=1,
        path_from=0, path_length=100, path_skew=0.5,
        path_nonlinearProb=0.1, path_sigmaFac=0.8,
        seed=0
    """,
    body="""
        group_prob <- as.numeric(group_prob)
        path_from <- as.numeric(path_from)
        path_length <- as.numeric(path_length)
        path_skew <- as.numeric(path_skew)
        sim <- splatSimulate(
            method=method, batchCells=batchCells, nGenes=nGenes,
            batch.facLoc=batch_facLoc, batch.facScale=batch_facScale,
            mean.rate=mean_rate, mean.shape=mean_shape,
            lib.loc=lib_loc, lib.scale=lib_scale, lib.norm=lib_norm,
            out.prob=out_prob,
            out.facLoc=out_facLoc, out.facScale=out_facScale,
            de.prob=de_prob, de.downProb=de_downProb,
            de.facLoc=de_facLoc, de.facScale=de_facScale,
            bcv.common=bcv_common, bcv.df=bcv_df,
            dropout.type=dropout_type, dropout.mid=dropout_mid,
            dropout.shape=dropout_shape,
            group.prob=group_prob,
            path.from=path_from, path.length=path_length, path.skew=path_skew,
            path.nonlinearProb=path_nonlinearProb, path.sigmaFac=path_sigmaFac,
            seed=seed
        )
        data <- sim %>%
            counts() %>%
            t()
        list(data=data, time=sim$Step, branch=sim$Group)
    """)


def SplatSimulate(
        method="paths",
        nBatches=1, batchCells=100,
        nGenes=10000,
        batch_facLoc=0.1, batch_facScale=0.1,
        mean_rate=0.3, mean_shape=0.6,
        lib_loc=11, lib_scale=0.2, lib_norm=False,
        out_prob=0.05,
        out_facLoc=4, out_facScale=0.5,
        de_prob=0.1, de_downProb=0.1,
        de_facLoc=0.1, de_facScale=0.4,
        bcv_common=0.1, bcv_df=60,
        dropout_type='none', dropout_prob=0.5,
        dropout_mid=0, dropout_shape=-1,
        group_prob=1,
        path_from=0, path_length=100, path_skew=0.5,
        path_nonlinearProb=0.1, path_sigmaFac=0.8,
        seed=None):
    """
    Global parameters
        nGenes - The number of genes to simulate.
        nCells - The number of cells to simulate.
        seed - Seed to use for generating random numbers.
    Batch parameters
        nBatches - The number of batches to simulate.
        batchCells - The number of cells in each batch.
        batch.facLoc - Location (meanlog) parameter for the batch effects factor log-normal distribution.
        batch.facScale - Scale (sdlog) parameter for the batch effects factor log-normal distribution.
    Mean parameters
        mean.shape - Shape parameter for the mean gamma distribution.
        mean.rate - Rate parameter for the mean gamma distribution.
    Library size parameters
        lib.loc - Location (meanlog) parameter for the library size log-normal distribution, or mean for the normal distribution.
        lib.scale - Scale (sdlog) parameter for the library size log-normal distribution, or sd for the normal distribution.
        lib.norm - Whether to use a normal distribution instead of the usual log-normal distribution.
    Expression outlier parameters
        out.prob - Probability that a gene is an expression outlier.
        out.facLoc - Location (meanlog) parameter for the expression outlier factor log-normal distribution.
        out.facScale - Scale (sdlog) parameter for the expression outlier factor log-normal distribution.
    Group parameters
        nGroups - The number of groups or paths to simulate.
        group.prob - The probabilities that cells come from particular groups.
    Differential expression parameters
        de.prob - Probability that a gene is differentially expressed in each group or path.
        de.loProb - Probability that a differentially expressed gene is down-regulated.
        de.facLoc - Location (meanlog) parameter for the differential expression factor log-normal distribution.
        de.facScale - Scale (sdlog) parameter for the differential expression factor log-normal distribution.
    Biological Coefficient of Variation parameters
        bcv.common - Underlying common dispersion across all genes.
        bcv.df - Degrees of Freedom for the BCV inverse chi-squared distribution.
    Dropout parameters
        dropout.type - Type of dropout to simulate.
        dropout.mid - Midpoint parameter for the dropout logistic function.
        dropout.shape - Shape parameter for the dropout logistic function.
    Differentiation path parameters
        path.from - Vector giving the originating point of each path.
        path.length - Vector giving the number of steps to simulate along each path.
        path.skew - Vector giving the skew of each path.
        path.nonlinearProb - Probability that a gene changes expression in a non-linear way along the differentiation path.
        path.sigmaFac - Sigma factor for non-linear gene paths.
    """
    if seed is None:
        seed = np.random.randint(2**16 - 1)
    if dropout_type == 'binomial':
        dropout_type = "none"
    else:
        dropout_prob = None
    np.random.seed(seed)

    data = _SplatSimulate(
        method=method,
        nBatches=nBatches, batchCells=batchCells,
        nGenes=nGenes,
        batch_facLoc=batch_facLoc, batch_facScale=batch_facScale,
        mean_rate=mean_rate, mean_shape=mean_shape,
        lib_loc=lib_loc, lib_scale=lib_scale, lib_norm=lib_norm,
        out_prob=out_prob,
        out_facLoc=out_facLoc, out_facScale=out_facScale,
        de_prob=de_prob, de_downProb=de_downProb,
        de_facLoc=de_facLoc, de_facScale=de_facScale,
        bcv_common=bcv_common, bcv_df=bcv_df,
        dropout_type=dropout_type, dropout_mid=dropout_mid,
        dropout_shape=dropout_shape,
        group_prob=group_prob,
        path_from=path_from, path_length=path_length, path_skew=path_skew,
        path_nonlinearProb=path_nonlinearProb, path_sigmaFac=path_sigmaFac,
        seed=seed)
    if dropout_prob is not None:
        data['data'] = np.random.binomial(n=data['data'], p=1 - dropout_prob,
                                          size=data['data'].shape)
    return data
