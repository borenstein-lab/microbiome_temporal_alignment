from scipy.spatial.distance import cdist
from scipy.spatial.distance import braycurtis
from scipy.spatial.distance import cdist
from scipy.spatial.distance import braycurtis
import numpy as np
import pandas as pd
from skbio import TreeNode
from skbio.diversity import beta_diversity






def pairwise_bray_curtis(sample1, sample2):
    return braycurtis(sample1,sample2)



def calculate_bray_curtis_table(OTU_table1, OTU_table2):
    matrix = cdist(OTU_table1.T, OTU_table2.T, metric='braycurtis')
    matrix = pd.DataFrame(matrix)
    matrix.columns = OTU_table2.columns
    matrix.index = OTU_table1.columns
    return matrix

def calculate_weighted_unifrac_table_wrapper(subject1,subject2,data):
    samples1 = data.counts_by_subject(subject1)
    samples2 = data.counts_by_subject(subject2)

    return calculate_weighted_unifrac_table(samples1,samples2,data)


def calculate_weighted_unifrac_table(samples1,samples2,data):
    ids1 = samples1.columns
    ids2 = samples2.columns
    both = pd.concat([samples1, samples2], axis=1)
    wu_dm = beta_diversity("weighted_unifrac", np.asarray(both.T), both.T.index, tree=data.tree, otu_ids=both.T.columns)
    unifrac = pd.DataFrame(wu_dm._data)
    unifrac.index = unifrac.columns = both.columns

    return unifrac.loc[ids1, ids2]


def calculate_time_delta_matrix(time_series1,time_series2, normalization="linear" ,threshold = 300, base=300):

    matrix = cdist(np.asarray(time_series1).reshape(-1, 1), np.asarray(time_series2).reshape(-1, 1))
    matrix = pd.DataFrame(matrix)
    matrix.columns = time_series2.index
    matrix.index = time_series1.index
    if normalization=="log":
        matrix = np.log(matrix)
        matrix=matrix/(np.log(base))
        matrix = matrix.clip(0,1)
    if normalization=="regular":
        matrix = matrix - matrix.min().min()
        matrix=matrix/matrix.max().max()

    if normalization == "linear":
        matrix = matrix / threshold
        matrix = matrix.clip(0,1)

    return matrix

class IntegratedDistanceMatrix(object):

    def __init__(self, w_t, w_d, T, D,standartize_matrices=False,normalize_matrices=False):
        """
            n - number of samples from first subject, m - number of samples from second subject
            :param
            T: n*m matrix, when T(i,j):= age of i'th sample of first subject -age of j'th sample of second subject
            D: a distance matrix, where D(i,j) is the distance between i'th sample of first subject and j'th sample of second subject,
            according to a given metric (e.g bray-curtis).
            subject: string, a subject number/ name. will be used to generate names for interpolated samples.
            :return:
            interpolated_taxa_table - a table with interpolated samples in columns and taxa (OTU, Genus etc.) on rows.
            interpolatd_ages - a list of the ages of interpolated samples, in the same order as in interpolated_taxa_table.
            """

        assert T.shape == D.shape

        self.w_t = w_t
        self.w_d = w_d
        if standartize_matrices:
            self.T = IntegratedDistanceMatrix.standartize(T)
            self.D = IntegratedDistanceMatrix.standartize(D)

        else:
            self.D = D
            self.T = T

        if normalize_matrices:
            self.T = IntegratedDistanceMatrix.normalize(self.T)
            self.D = IntegratedDistanceMatrix.normalize(self.D)
        T,D = pd.DataFrame(T), pd.DataFrame(D)
        matrix = (w_t*self.T).add(w_d*self.D)
        # print(w_t*self.T.shape,w_d*self.D.shape,"intermediate matrix shape:",matrix.shape)
        # print(matrix.shape)

        # Do I need to do this shift???
        # matrix-=matrix.min().min()
        # matrix/=matrix.max().max()
        self.IM=matrix

    @staticmethod
    def standartize(m):
        m -= np.mean(m)
        m /= np.std(m)
        return m

    @staticmethod
    def normalize(m):
        m-=m.min().min()
        m/=m.max().max()
        return m


