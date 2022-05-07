import numpy as np
import math
import pandas as pd
from matplotlib import pyplot as plt
import data_tables
import Config
import seaborn as sns


def gaussian_interpolation(window_size, taxa_table, subject, original_ages, normalize_time="no", day_interval=0, n=0, k_day_bins = False, k=30, plot=False):
    """
    :param
    taxa_table: a table with samples in columns and taxa (OTU, Genus etc.) on rows.
    original_ages: a list of the ages of samples in input, in the same order as the input in taxa_table
    subject: string, a subject number/ name. will be used to generate names for interpolated samples.
    :return:
    interpolated_taxa_table - a table with interpolated samples in columns and taxa (OTU, Genus etc.) on rows.
    interpolatd_ages - a list of the ages of interpolated samples, in the same order as in interpolated_taxa_table.
    """

    # taxa_table = samples_by_subject[subject]
    # original_ages = [int(ages_by_samples[sample][0]) for sample in taxa_table.columns]

    m = len(original_ages)
    # print(m)
    max_age = max(original_ages)
    min_age = min(original_ages)
    if (not n) and not k_day_bins:
        n = int((max_age - min_age) / day_interval)
    if k_day_bins:
        min_bin = math.ceil(min_age / k) * k
        max_bin = math.floor(max_age / k) * k
        n = int((max_bin-min_bin)/k)
    # print(max_age)
    weight_matrix = np.zeros(shape=(m, n), dtype=float)
    if normalize_time == "yes":
        # implement time normalization
        pass
    interpolated_ages = []
    tags = []
    for i in range(n):
        if(k_day_bins):
            inter_age = int(min_bin+i*k)
        else:
            inter_age = int(min_age + (i / (n - 1)) * (max_age - min_age))
        interpolated_ages.append(inter_age)
        tag = "I_{0}_{1}".format(inter_age, subject)
        tags.append(tag)
        # ages_by_samples[tag] = [str(inter_age)]
        # print(inter_age)
        for j in range(m):
            ref_age = original_ages[j]
            # print(inter_age,ref_age)
            dist = inter_age - ref_age
            weight_matrix[j][i] = math.exp(-(dist ** 2) / (window_size ** 2))
    # print(interpolated_ages)
    # print(weight_matrix)
    # print(taxa_table.shape)
    weights = pd.DataFrame(weight_matrix)
    weights = weights[~weights.isin([np.nan, np.inf, -np.inf]).any(1)]
    weight_matrix = weight_matrix / weight_matrix.sum(axis=0)
    interpolated_values = pd.DataFrame(np.matmul(np.asarray(taxa_table), weight_matrix))
    interpolated_values.index = taxa_table.index
    interpolated_values.columns = tags
    # samples_by_subject["I{0}".format(subject)] = interpolated_values
    # ages_by_subject["I{0}".format(subject)] = interpolated_ages
    # print(genus[[column for column in genus.columns if column not in interpolated_values.columns]].columns,"***",interpolated_values.columns)
    # OTU = pd.concat(
    #     [OTU[[column for column in OTU.columns if column not in interpolated_values.columns]], interpolated_values], 1)
    if (plot):
        plt.rcParams["figure.figsize"] = (15, 3)
        fig = plt.figure()
        plt.title("Interpolation of {0} with window size {1}".format(subject, window_size))
        gaussian_mean = original_ages[int(len(original_ages) / 2)]
        t = np.arange(min_age, max_age, 0.01)
        s = [math.exp(-((gaussian_mean - ti) ** 2) / (window_size ** 2)) for ti in t]
        s = s/np.sum(s)
        s = s/np.max(s) #scale the max to 1 for plot
        plt.scatter(original_ages, [0 for age in original_ages], c='b')
        plt.scatter(interpolated_ages, np.random.uniform(-0.05, 0.05, len(interpolated_ages)), c='r', marker='x')
        plt.plot(t, s,c="g")
        plt.show()

    return_metadata = pd.DataFrame(interpolated_ages)
    return_metadata.columns = ["Age_at_Collection"]
    return_metadata["SampleID"] = tags
    return_metadata["Subject_ID"] = subject
    return_metadata = return_metadata[["Subject_ID","SampleID","Age_at_Collection"]]

    # print(return_metadata)



    return interpolated_values, return_metadata


def interpolate_dataset(original_data, config,output_config=None):
    subjects = original_data.ages["Subject_ID"].unique()
    subject = subjects[0]
    taxa_table = original_data.abundance_by_subject(subject)
    ages = (original_data.ages_by_subject(subject))
    new_taxa, new_ages = gaussian_interpolation(config.window_size, taxa_table, subject, list(ages),
                                                day_interval=config.day_interval, n=config.number_of_days,
                                                k_day_bins=config.k_day_bins,k=config.k)
    for subject in subjects[1:]:
        # print(subject)
        taxa_table = original_data.abundance_by_subject(subject)
        ages = (original_data.ages_by_subject(subject))
        i_taxa, i_ages = gaussian_interpolation(config.window_size, taxa_table, subject, list(ages),
                                                day_interval=config.day_interval, n=config.number_of_days,
                                                k_day_bins=config.k_day_bins,k=config.k)
        new_taxa = pd.concat([new_taxa, i_taxa], axis=1)
        new_ages = pd.concat([new_ages, i_ages], axis=0)
    if type(output_config)==type(None):
        output_config = Config.Config(abundance_threshold=0, max_age=10000, min_max_age=0, window_size=10,
                                      top_predictors_for_age_analysis=0)

    return data_tables.ta_data(new_ages, abundance=new_taxa,config=output_config)