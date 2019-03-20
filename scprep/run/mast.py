import pandas as pd
from .r_function import RFunction


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
