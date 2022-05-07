import numpy as np
import math
import pandas as pd
from matplotlib import pyplot as plt
import random
import interpolation



def add_noise(taxa_table, ages, subject, window_size, mean = 0):
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
    shape = taxa_table.shape
    noise = np.random.normal(loc=mean, scale=window_size, size=shape)
    simulated_data = taxa_table + noise
    simulated_data.clip(lower=0,inplace=True)
    # print("simulated data shape:",simulated_data.shape)
    tags = ["S{0}{1}".format(age, subject) for age in ages]
    simulated_data.index = taxa_table.index
    simulated_data.columns = tags

    return_metadata = pd.DataFrame(ages)
    return_metadata.index = tags

    original_metadata = pd.DataFrame(ages)
    original_metadata.index = taxa_table.columns

    # return simulated_data, return_metadata, taxa_table, original_metadata
    return simulated_data, return_metadata

def shuffle_data(taxa_table, ages, subject, power):


    n = len(ages)
    indices = [i for i in range(n)]
    amount_to_shuffle = int(power*n)
    for i in range(amount_to_shuffle):
        index1 = np.random.randint(n)
        index2 = np.random.randint(n)
        indices[index1],indices[index2] = indices[index2],indices[index1]
    # print(indices)
    simulated_data = taxa_table.iloc[:,indices]
    # print("simulated data shape:",simulated_data.shape)
    shuffled_ages=[ages[indices[i]] for i in range(n)]
    tags = ["S{0}{1}".format(age, subject) for age in ages]
    simulated_data.index = taxa_table.index
    simulated_data.columns = tags

    # return_metadata = pd.DataFrame(shuffled_ages)
    return_metadata = pd.DataFrame(ages)
    return_metadata.index = tags

    original_metadata = pd.DataFrame(ages)
    original_metadata.index = taxa_table.columns
    # print(simulated_data.head(1))
    # print(return_metadata.head(1))

    # return simulated_data, return_metadata, taxa_table, original_metadata
    return simulated_data, return_metadata

def sub_sample_interpolation(taxa_table, ages, subject, amount, age_tags="first_n", interpolate=False, interpolation_window=8):
    print(ages)
    n = len(ages)
    print(n)
    indices = [i for i in range(n)]
    sample = sorted(random.sample(indices, amount))
    print(sample)
    simulated_data_seed = taxa_table.iloc[:, sample]
    if age_tags=="actual":
        simulated_ages_seed = [ages[sample[i]] for i in range(len(sample))]
    if age_tags=="first_n":
        simulated_ages_seed = [ages[i] for i in range(len(sample))]
    else:
        simulated_ages_seed = age_tags
    # print(simulated_ages_seed)
    tags = ["S{0}{1}".format(age, subject) for age in simulated_ages_seed]
    simulated_data_seed.index = taxa_table.index
    simulated_data_seed.columns = tags
    if interpolate:
        simulated_data, return_metadata = interpolation.gaussian_interpolation(interpolation_window, simulated_data_seed, subject, simulated_ages_seed, normalize_time="no", n=n)
    else:
        return_metadata = pd.DataFrame(simulated_ages_seed)
        return_metadata.index = tags
        simulated_data = simulated_data_seed
    # return_metadata = pd.DataFrame(simulated_ages_seed)
    # return_metadata.index = tags

    original_metadata = pd.DataFrame(ages)
    original_metadata.index = taxa_table.columns
    # print(simulated_data.head(1))
    # print(return_metadata.head(1))

    return simulated_data, return_metadata, taxa_table, original_metadata


def add_tail(taxa_table, ages, subject, length, day_interval, data):
    n = len(ages)
    max_age = max(ages)
    tail_ages = [max_age+i*day_interval for i in range(1,length+1)]
    indices = [i for i in range(data.counts.shape[1])]
    sample = sorted(random.sample(indices, length))
    print(sample)
    tail = data.counts.iloc[:, sample]
    simulated_data_seed = pd.concat([taxa_table,tail],axis=1)
    print(taxa_table.shape)
    print(simulated_data_seed.shape)
    print(len(ages), len(tail_ages))
    all_ages = ages + tail_ages
    tags = ["S{0}{1}".format(age, subject) for age in all_ages]
    simulated_data_seed.index = taxa_table.index
    simulated_data_seed.columns = tags
    return_metadata = pd.DataFrame(all_ages)
    return_metadata.index = tags
    simulated_data = simulated_data_seed
    # return_metadata = pd.DataFrame(simulated_ages_seed)
    # return_metadata.index = tags

    original_metadata = pd.DataFrame(ages)
    original_metadata.index = taxa_table.columns
    # print(simulated_data.head(1))
    # print(return_metadata.head(1))

    return simulated_data, return_metadata, taxa_table, original_metadata