import pandas as pd
import numpy as np
import random
import distance
import alignment
from matplotlib import pyplot as plt
import seaborn as sns

def crop_by_ages(subject1,subject2,data,match_length=False):
    samples1, samples2 = data.sample_ids_by_subject(subject1), data.sample_ids_by_subject(subject2)
    ages1, ages2 = data.get_ages_of_samples(samples1), data.get_ages_of_samples(samples2)
    min_age = max(min(ages1),min(ages2))
    max_age = min(max(ages1),max(ages2))
    ages2 = ages2[(ages2>=min_age) & (ages2<=max_age)]
    ages1 = ages1[(ages1>=min_age) & (ages1<=max_age)]
    if match_length:
        len_to_crop = min(len(ages1),len(ages2))
        ages1,ages2 = ages1[:len_to_crop],ages2[:len_to_crop]
#     print(ages1)
#     print(ages2)
#     print("len:",len(ages1),len(ages2))
    samples1 = list(ages1.index)
    samples2 = list(ages2.index)
    return samples1, ages1, samples2, ages2



def calc_pairwise_alignment_matrix(data, distance_matrix, plot=False, w_t=0, w_d=1, sub_sample=0,step_pattern="symmetric2",window_type="none",window_args={},crop=False):
    subjects = data.all_subjects()
    n = len(subjects)
    results = pd.DataFrame(np.zeros([n, n]))
    results.index = results.columns = subjects
    for i, subject1 in enumerate(subjects):
        for subject2 in subjects[i + 1:]:
            if crop:
                samples1, ages1, samples2, ages2 = crop_by_ages(subject1,subject2,data,match_length=False)
            else:
                samples1, samples2 = data.sample_ids_by_subject(subject1), data.sample_ids_by_subject(subject2)
            if sub_sample > 0:
                start_index1 = random.randint(0, len(samples1) - sub_sample)
                samples1 = samples1[start_index1:start_index1 + sub_sample]
                start_index2 = random.randint(0, len(samples2) - sub_sample)
                samples2 = samples2[start_index2:start_index2 + sub_sample]
            if not crop:
                ages1, ages2 = data.get_ages_of_samples(samples1), data.get_ages_of_samples(samples2)
            time_delta = distance.calculate_time_delta_matrix(ages1, ages2, normalization="linear", threshold=200)
            dist = distance_matrix.loc[samples1, samples2]
            integrated_distance = distance.IntegratedDistanceMatrix(w_t=w_t, w_d=w_d, T=time_delta, D=dist,
                                                                    standartize_matrices=False,
                                                                    normalize_matrices=False).IM
            alignment_object = alignment.align(x=None, y=None, distance_matrix=integrated_distance, how="global",
                                               threshold=0,step_pattern=step_pattern,window_type=window_type,window_args=window_args)
            if plot:
                fig, ax = plt.subplots(figsize=(3, 3))
                ii, jj = alignment_object.index1, alignment_object.index2
                sns.heatmap(dist)
                plt.plot(jj,ii)
                plt.title(f"{subject1} vs. {subject2}, distance={alignment_object.normalizedDistance}")
                plt.show()
            results.loc[subject1, subject2] = results.loc[subject2, subject1] = alignment_object.normalizedDistance
    return results

def calc_alternative_distance_matrix(data, distance_matrix, plot=False, w_t=0, sub_sample=0,kind="diagonal",crop=False):
    subjects = data.all_subjects()
    n = len(subjects)
    results = pd.DataFrame(np.zeros([n, n]))
    results.index = results.columns = subjects
    for i, subject1 in enumerate(subjects):
        for subject2 in subjects[i + 1:]:
            if crop:
                samples1, ages1, samples2, ages2 = crop_by_ages(subject1,subject2,data,match_length=False)
            else:
                samples1, samples2 = data.sample_ids_by_subject(subject1), data.sample_ids_by_subject(subject2)
            if sub_sample > 0:
                start_index1 = random.randint(0, len(samples1) - sub_sample)
                samples1 = samples1[start_index1:start_index1 + sub_sample]
                start_index2 = random.randint(0, len(samples2) - sub_sample)
                samples2 = samples2[start_index2:start_index2 + sub_sample]
            if not crop:
                ages1, ages2 = data.get_ages_of_samples(samples1), data.get_ages_of_samples(samples2)
            time_delta = distance.calculate_time_delta_matrix(ages1, ages2, normalization="linear", threshold=200)
            dist = distance_matrix.loc[samples1, samples2]
            integrated_distance = distance.IntegratedDistanceMatrix(w_t=w_t, w_d=1, T=time_delta, D=dist,
                                                                    standartize_matrices=False,
                                                                    normalize_matrices=False).IM
            if kind=="diagonal":
                diagonal_distance = np.mean(np.diag(integrated_distance))
                results.loc[subject1, subject2] = results.loc[subject2, subject1] = diagonal_distance
            if kind=="min":
                min_distance = integrated_distance.min().min()
                results.loc[subject1, subject2] = results.loc[subject2, subject1] = min_distance
            if kind=="mean":
                mean_distance = integrated_distance.mean().mean()
                results.loc[subject1, subject2] = results.loc[subject2, subject1] = mean_distance


    return results


# twin_pairs =
def get_family(full_metadata, subject):
    return full_metadata[full_metadata["PersonID"] == subject]["FamilyID"].iloc[0]


def get_subject(subject):
    return subject[:-2]


def get_case_control(metadata, subject):
    return metadata[metadata["Subject_ID"] == subject]["Case_Control"].iloc[0]


def calc_twin_no_twin(metadata, alignment_distances,split_data=True):
    twin_distances = []
    non_twin_distances = []
    same_subject_distances = []
    subjects = alignment_distances.columns
    for i, subject1 in enumerate(subjects):
        for subject2 in subjects[i + 1:]:
            if subject1 != subject2:
                dist = alignment_distances.loc[subject1, subject2]
                if split_data:
                    if get_subject(subject1) == get_subject(subject2):
                        same_subject_distances += [dist]
                    elif get_family(metadata, subject1[:-2]) == get_family(metadata, subject2[:-2]):
                            twin_distances += [dist]
                    else:
                        non_twin_distances += [dist]
                else:
                    if get_family(metadata, subject1) == get_family(metadata, subject2):
                        twin_distances += [dist]
                    else:
                        non_twin_distances += [dist]
    return same_subject_distances, twin_distances, non_twin_distances


def calc_same_subject_distances(alignment_distances):
    different_subject_distances = []
    same_subject_distances = []
    subjects = alignment_distances.columns
    for i, subject1 in enumerate(subjects):
        for subject2 in subjects[i + 1:]:
            if subject1 != subject2:
                dist = alignment_distances.loc[subject1, subject2]
                if subject1[:-1]==subject2[:-1]:
                    same_subject_distances += [dist]
                else:
                    different_subject_distances += [dist]
    return same_subject_distances, different_subject_distances


def calc_case_control(metadata, alignment_distances):
    case_distances = []
    control_distances = []
    between_distances = []
    subjects = metadata["Subject_ID"].unique()
    for i, subject1 in enumerate(subjects):
        for subject2 in subjects[i + 1:]:
            dist = alignment_distances.loc[subject1, subject2]
            if get_case_control(metadata, subject1) == get_case_control(metadata, subject2) == "case":
                case_distances += [dist]
            if get_case_control(metadata, subject1) == get_case_control(metadata, subject2) == "control":
                control_distances += [dist]
            if get_case_control(metadata, subject1) != get_case_control(metadata, subject2):
                between_distances += [dist]
    return case_distances, control_distances, between_distances



