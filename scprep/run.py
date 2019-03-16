try:
    from rpy2.robjects.packages import STAP
    import rpy2.robjects.numpy2ri
    import rpy2.robjects  # NOQA
except ImportError:
    pass
import numpy as np
import pandas as pd
from decorator import decorator


@decorator
def _with_rpy2(fun, *args, **kwargs):
    try:
        rpy2
    except NameError:
        raise ImportError(
            "rpy2 not found. "
            "Please install it with e.g. `pip install --user rpy2`")
    return fun(*args, **kwargs)


class RFunction(object):

    @_with_rpy2
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
            self._function = getattr(STAP(function_text, self.name), self.name)
            rpy2.robjects.numpy2ri.activate()
            return self._function

    def is_r_object(self, obj):
        return "rpy2.robjects" in str(type(obj))

    def convert(self, robject):
        if self.is_r_object(robject):
            if isinstance(robject, rpy2.robjects.vectors.ListVector):
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
                robject = rpy2.robjects.numpy2ri.ri2py(robject)
                if self.is_r_object(robject):
                    # try regular conversion
                    robject = rpy2.robjects.conversion.ri2py(robject)
        return robject

    def __call__(self, *args, **kwargs):
        robject = self.function(*args, **kwargs)
        return self.convert(robject)

# TODO: Figure out how to check if the nescessary package is already installed
# TODO: Make rpy2 not a depedency
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
