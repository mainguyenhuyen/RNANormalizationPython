# Programming project 

# Import required packages
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri

# Convert pandas.DataFrames to R dataframes automatically.
pandas2ri.activate()

class RNARawCount:
    '''
    Description: This is a Python class specialized in RNA sequencing (RNA-seq) normalization.  

    '''
    
    def __init__(self, fileLocation, lowCountLimit = 10):
        '''
        The initalize funtion.
        
        Parameters
        ----------
        fileLocation : str
            DESCRIPTION. the directory where the count matrix is stored. 
        lowCountLimit : int, optional
            DESCRIPTION. the number as the lowest limit for low count filtration. The default value is 10.

        Raises
        ------
        ValueError
            DESCRIPTION. Two situations can raise ValueError: input file is not a comma-separated values or tab-separated or text; or input file does not have gene length column with the name 'width'.

        Returns
        -------
        None.

        '''
        # Determine the file extension
        if fileLocation.endswith('.csv'):
            fileType = 'csv'
        elif fileLocation.endswith('.txt'):
            fileType = 'txt'
        elif fileLocation.endswith('.tsv'):
            fileType = 'tsv'
        else:
            raise ValueError('Unsupported file type')
        # Read the file into a data frame
        if fileType == 'csv':
            dataCount = pd.read_csv(fileLocation, sep=',', index_col = 0)
        elif fileType == 'txt' or 'tsv':
            dataCount = pd.read_csv(fileLocation, sep='\t', index_col = 0)
        if "width" not in dataCount.columns:
            raise ValueError('File was unsuccessfully read due to incorrect format')
        # Pre-processing
        # Change the data type to np.float64 for calculation later
        dataCount = dataCount.astype(np.float64)
        # Create a data frame contains only reads count, no gene length 
        readsOnly = dataCount.drop('width', axis = 1)
        # Filter low count genes by taking genes with total count of all samples more than lowCountLimit
        self.dataCount = dataCount[readsOnly.sum(axis = 1) > lowCountLimit]
        # Re-create a data frame contains only reads count, no gene length after the filtering step 
        self.readsOnly = self.dataCount.drop('width', axis = 1)
        # Extract the gene length of the count data
        self.geneLength = self.dataCount.pop('width') 
        # Compute the total number of read mapped over all genes per sample
        self.totalMappedReads = self.readsOnly.sum(axis=0)
    def __len__(self):
        '''
        Returning the number of rows in the count matrix.

        Returns:
        -------
        int
            DESCRIPTION. the number of rows or genes in the count matrix.

        '''
        return len(self.dataCount.index)
    def plotBeforeNormaliztion(self):
        '''
        Generating a histogram plot of the data from the count matrix.

        Returns
        -------
        None.

        '''
        plt.hist(self.readsOnly, bins = 50)
        plt.show()
    # Normalization
    # Choose which type of test is used
    def TPM(self):
        '''
        Computing normalized counts by TPM calculation.

        Returns
        -------
        TPM_result : dataframe
            DESCRIPTION. A data frame containing the normalized counts by TPM calculation.

        '''
        # Compute reads multiply with 1e6
        readsMultiply1e6 = self.readsOnly*1e6
        # Compute TPM
        firstDivisionTPM = readsMultiply1e6.divide(self.geneLength, axis = 0)
        RPK_TPM = self.readsOnly.divide(self.geneLength, axis = 0)
        totalRPK = RPK_TPM.sum(axis = 0)
        TPM_result = firstDivisionTPM.divide(totalRPK, axis = 1)
        return TPM_result
    def CPM(self):
        '''
        Computing normalized counts by CPM calculation.

        Returns
        -------
        CPM_result : dataframe
            DESCRIPTION. A data frame containing the normalized counts by CPM calculation.

        '''
        # Compute reads multiply with 1e9
        readsMultiply1e6 = self.readsOnly*1e6
        # Compute cpm
        CPM_result = readsMultiply1e6.divide(self.totalMappedReads, axis = 1)
        return CPM_result
    def RPKM(self):
        '''
        Computing normalized counts by RPKM calculation.

        Returns
        -------
        res : dataframe
            DESCRIPTION. A data frame containing the normalized counts by RPKM calculation.

        '''
        # Compute reads multiply with 1e6
        readsMultiply1e9 = self.readsOnly*1e9
        # Compute RPKM
        firstDivisionRes = readsMultiply1e9.divide(self.geneLength, axis = 0)
        res = firstDivisionRes.divide(self.totalMappedReads, axis = 1)
        return res
    def FPKM(self):
        '''
        Computing normalized counts by FPKM calculation.

        Returns
        -------
        res : dataframe
            DESCRIPTION. A data frame containing the normalized counts by FPKM calculation.

        '''
        # Compute reads multiply with 1e6
        readsMultiply1e9 = self.readsOnly*1e9
        # Compute RPKM
        firstDivisionRes = readsMultiply1e9.divide(self.geneLength, axis = 0)
        res = firstDivisionRes.divide(self.totalMappedReads, axis = 1)
        return res        
    def logarithmTransformation(self):
        '''
        Computing normalized counts by logarithm transformation of (x + 1).

        Returns
        -------
        logNormalization : dataframe
            DESCRIPTION. A data frame containing the normalized counts by log transformation.

        '''
        # Compute log normalization 
        logNormalization = np.log(self.readsOnly+1)
        return logNormalization
    def DESeq2MedianOfRatios(self):
        '''
        Computing normalized counts by DESeq2 median of ratios method.

        Returns
        -------
        normCounts : dataframe
            DESCRIPTION. A data frame containing the normalized counts by DESeq2 median of ratios method.

        '''
        # Remove genes if in any sample it has 0 count
        readsNoZero = self.readsOnly[(self.readsOnly != 0).all(axis=1)]
        # Extract the number of samples
        root = len(readsNoZero.columns)
        #compute row-wise geometric mean
        productOfCountPerGene = readsNoZero.prod(axis = 1)
        pseudoRef = productOfCountPerGene**(1/root)
        # Compute ratio of each sample
        ratioOfSamplePerRef = readsNoZero.divide(pseudoRef, axis = 0)
        # Compute normalization factor
        normFactor = np.median(ratioOfSamplePerRef, axis = 0)
        # Compute normalized counts
        normCounts = readsNoZero.divide(normFactor)
        return normCounts 
    def rlog(self, sampleType, design):
        '''
        Computing normalized counts by rlog function from DESeq2 R package.

        Parameters
        ----------
        sampleType : list
            DESCRIPTION. a list containing the sample types of all samples.
        design : str
            DESCRIPTION. a string containing the formula of sample types.

        Returns
        -------
        res : dataframe
            DESCRIPTION. A data frame containing the normalized counts by rlog function from DESeq2 R package.

        '''
        # Import DESeq2 and SummarizedExperiment library
        deseq = importr('DESeq2')
        sumdeseq = importr('SummarizedExperiment')
        # Create parameters for deseq matrix
        new = {}
        new['sampleType'] = robjects.StrVector(sampleType)
        colData = robjects.DataFrame(new)
        design = robjects.Formula(design)
        # Create a deseq count matrix input
        dds = deseq.DESeqDataSetFromMatrix(countData = self.readsOnly, colData = colData, design = design)
        normCounts = deseq.rlog(dds) 
        res = sumdeseq.assay(normCounts)
        res = pd.DataFrame(res, index = self.readsOnly.index.values, columns = self.readsOnly.columns)
        return res
    def vst(self, sampleType, design):
        '''
        Computing normalized counts by VST function from DESeq2 R package.

        Parameters
        ----------
        sampleType : list
            DESCRIPTION. a list containing the sample types of all samples
        design : str
            DESCRIPTION. a string containing the formula of sample types.

        Returns
        -------
        res : dataframe
            DESCRIPTION. A data frame containing the normalized counts by VST function from DESeq2 R package.

        '''
        # Import DESeq2 and SummarizedExperiment library
        deseq = importr('DESeq2')
        sumdeseq = importr('SummarizedExperiment')
        # Create parameters for deseq matrix
        new = {}
        new['sampleType'] = robjects.StrVector(sampleType)
        colData = robjects.DataFrame(new)
        design = robjects.Formula(design)
        # Create a deseq count matrix input
        dds = deseq.DESeqDataSetFromMatrix(countData = self.readsOnly, colData = colData, design = design)
        # Compute vst normalization
        normCounts = deseq.varianceStabilizingTransformation(dds)
        # Extract the normalized counts
        res = sumdeseq.assay(normCounts)
        res = pd.DataFrame(res, index = self.readsOnly.index.values, columns = self.readsOnly.columns)
        return res