import numpy as np
import pandas as pd
import warnings

from . import utils

rpy2 = utils.try_import("rpy2")
robjects = utils.try_import("rpy2.robjects")

_formatwarning = warnings.formatwarning


def _quiet_rwarning(message, category, *args, **kwargs):
    if category == rpy2.rinterface.RRuntimeWarning:
        return str(message) + '\n'
    else:
        return _formatwarning(message, category, *args, **kwargs)


class RFunction(object):

    @utils._with_pkg("rpy2")
    def __init__(self, name, args, setup, body, quiet_setup=True):
        self.name = name
        self.args = args
        self.setup = setup
        self.body = body
        if quiet_setup:
            self.setup = """
                suppressPackageStartupMessages(suppressMessages(
                    suppressWarnings({{
                        {setup}
                    }})))""".format(setup=self.setup)

    @property
    def function(self):
        try:
            return self._function
        except AttributeError:
            function_text = """
            {setup}
            {name} <- function({args}) {{
              {body}
            }}
            """.format(setup=self.setup, name=self.name,
                       args=self.args, body=self.body)
            self._function = getattr(robjects.packages.STAP(
                function_text, self.name), self.name)
            robjects.numpy2ri.activate()
            return self._function

    def is_r_object(self, obj):
        return "rpy2.robjects" in str(type(obj))

    def convert(self, robject):
        if self.is_r_object(robject):
            if isinstance(robject, robjects.vectors.ListVector):
                names = self.convert(robject.names)
                if names is rpy2.rinterface.NULL or \
                        len(names) != len(np.unique(names)):
                    # list
                    robject = [self.convert(obj) for obj in robject]
                else:
                    # dictionary
                    robject = {name: self.convert(
                        obj) for name, obj in zip(robject.names, robject)}
            else:
                # try numpy first
                robject = robjects.numpy2ri.ri2py(robject)
                if self.is_r_object(robject):
                    # try regular conversion
                    robject = robjects.conversion.ri2py(robject)
                if robject is rpy2.rinterface.NULL:
                    robject = None
        return robject

    def __call__(self, *args, **kwargs):
        # monkey patch warnings
        warnings.formatwarning = _quiet_rwarning
        robject = self.function(*args, **kwargs)
        robject = self.convert(robject)
        warnings.formatwarning = _formatwarning
        return robject


# TODO: Figure out how to check if the nescessary package is already installed
# TODO: Is keeping these scripts in this file the best way to maintain code?
# TODO: Make geneson optional for MAST, adding option for releveling the
# condition factor
_Monocle2 = RFunction(
    name="run_monocle",
    setup="library(monocle)",
    args="data",
    body="""
    data <- t(data)
    colnames(data) <- 1:ncol(data)
    rownames(data) <- 1:nrow(data)
    fd <- new("AnnotatedDataFrame", data = as.data.frame(rownames(data)))
    data <- newCellDataSet(data,phenoData=NULL,featureData = fd,
                         expressionFamily=uninormal())
    varLabels(data@featureData) <- 'gene_short_name'
    data <- estimateSizeFactors(data)
    data_reduced <- suppressMessages(suppressWarnings(reduceDimension(data,
                                   max_components=2,
                                   reduction_method='DDRTree',
                                   norm_method="none", scaling=FALSE)))

    # 2D embedding
    cell_embedding = reducedDimS(data_reduced)
    t(cell_embedding)""")

_MAST = RFunction(
    name="run_MAST",
    setup="""library(MAST)
             library(data.table)""",
    args="data, gene_names, condition_labels",
    body="""
    data <- t(data)

    fdat <- data.frame(primerid = factor(gene_names))

    sca <- FromMatrix(data, fData = fdat)

    # calculate number of genes detected per cell, called DetRate
    cdr2 <- colSums(assay(sca)>0)
    colData(sca)$cngeneson <- scale(cdr2, scale = FALSE) # zscore

    colData(sca)$condition <- factor(condition_labels)
    zlmCond <- zlm(~condition + cngeneson, sca)
    contr <- paste('condition', tail(levels(factor(condition_labels)), n=1), sep='')

    summaryCond <- summary(zlmCond, doLRT=contr)

    summaryDt <- summaryCond$datatable
    fcHurdle <- merge(summaryDt[contrast==contr & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                    summaryDt[contrast==contr & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients

    fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
    fcHurdleSig <- merge(fcHurdle, as.data.table(mcols(sca)), by='primerid')
    setorder(fcHurdleSig, fdr)
    fcHurdleSig <- t(fcHurdleSig)
    return(fcHurdleSig)""")

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
            dropout.type=dropout_type, dropout.mid=dropout_mid, dropout.shape=dropout_shape,
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


def MAST(data, gene_names, condition_labels):
    """Short function descriptor

    Takes a data matrix and performs pairwise differential expression analysis
    using a Hurdle model as implemented in [MAST](https://github.com/RGLab/MAST/).
    The current implementation uses the Cell Detection Rate (# non-zero genes per cell)
    as a factor in the analysis because this was found to be higher performing in
    a comparison of differential expression tools (https://www.ncbi.nlm.nih.gov/pubmed/29481549).

    Parameters
    ----------
    data : array-like, shape=[n_samples, n_features]
        Input data from both samples
    gene_names : array-like, shape=[n_features]
        Names of each gene - must be unique
    condition_labels : array-like, shape=[n_samples]
        A vector of condition labels associated with each sample.
        E.g. `['Expt', 'Expt', ... , 'Ctrl']`.
        Note that `len(np.unique(condition_labels))` must be `2`.

    Returns
    -------
    results : pandas.DataFrame
        A table of results with columns `Pr(>Chisq)`, `coef`, `ci.hi`, `ci.lo`, and	`fdr`.
        `Pr(>Chisq)`: the p-value associated with the gene
        `coef`: the estimated log-Fold-Change (logFC) associated with the Hurdle model
        `ci.hi`: upper bound of the confidence interval for logFC
        `ci.lo`: lower bound of the confidence interval for logFC
        `fdr`: false-discovery rate
        The reccomended significant threshold from this table is:
            fdr <= 0.05
            coef >= log2(1.5)
    Examples
    --------
    >>> import scprep
    >>> data = scprep.io.load_csv("my_data.csv")
    >>> data_ln = scprep.normalize.library_size_normalize(data)
    >>> data_log = scprep.transform.log(data, base=2)
    >>> cond = np.hstack([np.tile('cond1', ncells_in_cond1), np.tile('cond2', ncells_in_cond2)])
    >>> results = scprep.run.run_MAST(data, gene_names=data.columns, condition=cond)
    """
    results = pd.DataFrame.from_records(_MAST(
        data, gene_names=gene_names, condition_labels=condition_labels), index='primerid')
    results.index.names = ['gene_name']
    return results


def Monocle2(data):
    """Short function descriptor

    Long function descriptor

    Parameters
    ----------
    data : array-like, shape=[n_samples, n_features]
        Input data from both samples

    Returns
    -------
    embedding : data type
        Description

    Examples
    --------
    >>> import scprep
    >>> data = scprep.io.load_csv("my_data.csv")
    >>> data_ln = scprep.normalize.library_size_normalize(data)
    >>> data_log = scprep.transform.log(data)
    >>> results = scprep.run.Monocle(data_log)
    """
    return _Monocle2(data)


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
