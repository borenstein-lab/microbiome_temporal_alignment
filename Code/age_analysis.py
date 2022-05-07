import sys

import matplotlib.cm as cm
import scipy
import sklearn
import seaborn as sns
import alignment
import pandas as pd
import numpy as np
import random
from matplotlib import pyplot as plt
from sklearn.ensemble import RandomForestRegressor
import time



def predict_ages_single_pair(ages1, ages2, ii, jj):
    df1, df2 = pd.DataFrame(ages1[ii]).reset_index(drop=True), pd.DataFrame(ages2[jj]).reset_index(drop=True)
    merged = pd.concat([df1, df2], axis=1)
    merged.index = ages1[ii].index
    merged = merged.groupby(merged.index).mean()
    merged = merged.reindex(ages1.index)
    merged.columns = ["drop", "predicted_age"]
    merged = merged.drop("drop", axis=1)
    return merged


def predict_ages(samples1, samples2, data, distance_matrix, measurements, config, how="all", alignemnt_method="open_begin_end", input_type="distance"):
    ages1, ages2 = data.get_ages_of_samples(samples1), data.get_ages_of_samples(samples2)
    if input_type=="distance":
        dist = distance_matrix.loc[samples1, samples2]
        dist.index = ages1
        dist.columns = ages2
        if config.align_threshold_quantile:
            thresh = np.mean(dist.quantile(q=config.align_threshold_quantile))
        else:
            thresh = config.align_threshold_value
        result = alignment.align(x=None,y=None, distance_matrix=dist, how=alignemnt_method, threshold=thresh, step_pattern="symmetric2",input_type="distance")
    elif input_type=="measurements":
        measurements1 = measurements[samples1]
        measurements2 = measurements[samples2]
        result = alignment.align(x=measurements1, y=measurements2, distance_matrix=None, how=alignemnt_method, threshold=0, step_pattern="symmetric2",input_type="measurements")
    # print(result.dtw_object.normalizedDistance)
    ii, jj = result.ii, result.jj
    df1, df2 = pd.DataFrame(ages1[ii]).reset_index(drop=True), pd.DataFrame(ages2[jj]).reset_index(drop=True)
    merged = pd.concat([df1, df2], axis=1)
    merged.index = ages1[ii].index
    merged = merged.groupby(merged.index).mean()
    merged = merged.reindex(ages1.index)
    merged.columns = ["drop", "predicted_age"]
    merged["alignment_score"]=result.dtw_object.normalizedDistance
    merged = merged.drop("drop", axis=1)

    return merged


def plot_alternative_results(results,compact=True,title="placeholder"):
    # results = results[["sample_ids","true","rf_predictions","1","2","3","4","5"]]
    for_plot = [results[["true",column]] for column in results.columns[3:]]
    if compact:
        ax = plot_age_prediction_compact(for_plot, label="DI_filtered_results", color="g",title=title)
    else:
        ax = plot_age_prediction(for_plot,label="DI_filtered_results",color="g")
    # add_rf_to_plot(ax,pd.DataFrame(results[["true","rf_predictions"]]))
    plt.show()

# max_start_index must be >= n!
def predict_ages_for_n_contigs(n, max_start_index, iterations, data, config, distance_matrix, measurements, input_type, how="all"):
    true = []
    predicted = []
    colors = cm.rainbow(np.linspace(0, 1, iterations))
    subjects = list(data.ages["Subject_ID"].unique())
    for i in range(iterations):
        # all_predictions = []
        all_predictions = pd.DataFrame(columns=["predicted_age", "alignment_score"], dtype=float)
        subject1 = random.choice(subjects)
        all_samples = data.sample_ids_by_subject(subject1)
        #         start_index = random.randint(0,len(all_samples)-n)
        start_index = random.randint(0, len(all_samples) - max_start_index)
        samples1 = all_samples[start_index:start_index + n]
        ages1 = data.get_ages_of_samples(samples1)
        for subject2 in [subject for subject in subjects if subject != subject1]:
            samples2 = data.sample_ids_by_subject(subject2)
            predictions = predict_ages(samples1, samples2, data, distance_matrix, measurements, input_type,config)
            # predictions = list(predict_ages(samples1, samples2, data, distance_matrix, config).iloc[:, 0])
            all_predictions = pd.concat([all_predictions,predictions])
        all_predictions.index.rename("index", inplace=True)
        if config.top_predictors_for_age_analysis != 0:
            print(all_predictions)
            print(all_predictions_filtered)
            threshold_score = sorted(set(all_predictions["alignment_score"]))[config.top_predictors_for_age_analysis]
            all_predictions_filtered = all_predictions[all_predictions["alignment_score"] < threshold_score]
            all_predictions = all_predictions_filtered
        # print(all_predictions.drop("alignment_score", axis=1).groupby("index").mean())
        mean_predictions = list(all_predictions.drop("alignment_score", axis=1).groupby("index").mean().iloc[:,0])
        # print(mean_predictions)
        if how == "first":
            true += [list(ages1)[0]]
            predicted += [mean_predictions[0]]
        if how == "all":
            true += list(ages1)
            predicted += mean_predictions
    #         print(true,predicted)
    #         plt.scatter(list(ages1), mean_predictions)
    results = pd.DataFrame([true, predicted]).T
    results = results.dropna()
    results.index = [i for i in range(len(results[0]))]
    # print(results)

    return results


def alternative_predict_ages_for_n_contigs(n_list, max_start_index, iterations, data, config, distance_matrix, measurements, input_type, how="all"):
    true = []
    predicted = [[] for n in range(len(n_list)+1)]
    rf_predictions = []
    sample_ids = []
    subjects = list(data.ages["Subject_ID"].unique())
    start_time = time.perf_counter()
    for i in range(iterations):
        print(f"i={i}, time from start={time.perf_counter()-start_time}", end='\r')
        sys.stdout.flush()
        # all_predictions = []
        # all_predictions = pd.DataFrame(columns=["predicted_age", "alignment_score"], dtype=float)
        subject1 = random.choice(subjects)
        all_samples = data.sample_ids_by_subject(subject1)
        #         start_index = random.randint(0,len(all_samples)-n)
        start_index = random.randint(0, len(all_samples) - max_start_index)
        start_sample = all_samples[start_index]
        sample_ids += [start_sample]
        true_age = data.age_of_sample(start_sample).iloc[0]
        # print(true_age)
        if how == "first":
            true += [true_age]
        for n in n_list:
            all_predictions = pd.DataFrame(columns=["predicted_age", "alignment_score"], dtype=float)
            samples1 = all_samples[start_index:start_index + n]
            # print(n,samples1)
            ages1 = data.get_ages_of_samples(samples1)
            for subject2 in [subject for subject in subjects if subject != subject1]:
                samples2 = data.sample_ids_by_subject(subject2)
                predictions = predict_ages(samples1, samples2, data, distance_matrix=distance_matrix, config=config, measurements=measurements, input_type=input_type)
                if how == "first":
                    # print(samples1[0])
                    predictions = predictions[predictions.index==start_sample]
                # predictions = list(predict_ages(samples1, samples2, data, distance_matrix, config).iloc[:, 0])
                all_predictions = pd.concat([all_predictions,predictions])
            all_predictions.index.rename("index", inplace=True)
            if config.top_predictors_for_age_analysis != 0:
                threshold_score = sorted(set(all_predictions["alignment_score"]))[config.top_predictors_for_age_analysis]
                all_predictions_filtered = all_predictions[all_predictions["alignment_score"] < threshold_score]
                all_predictions = all_predictions_filtered
            # print(all_predictions.drop("alignment_score", axis=1).groupby("index").mean())
            mean_predictions = list(all_predictions.drop("alignment_score", axis=1).groupby("index").mean().iloc[:,0])
            # print(mean_predictions)
            if how == "first":
                predicted[n] += [mean_predictions[0]]
                # print(predicted)
            if how == "all":
                predicted[n] += mean_predictions
        #         print(true,predicted)
        #         plt.scatter(list(ages1), mean_predictions)

        regr = predict_ages_with_rf(data, all_samples)[0]
        rf_prediction = regr.predict(data.abundance[start_sample].values.reshape(1, -1))
        rf_predictions += list(rf_prediction)

    results_dict = {str(n): predicted[n] for n in n_list}
    results_dict["sample_ids"] = sample_ids
    results_dict["true"]=true
    results_dict["rf_predictions"]=rf_predictions
    results = pd.DataFrame(results_dict)
    results = results.dropna()
    results = results[["sample_ids","true","rf_predictions"]+[str(n) for n in n_list]]
    results.set_index("sample_ids")
    # results.index = [i for i in range(len(results[0]))]
    # print(results)

    return results


def run_predictions_on_contig_len_list(contig_lengths ,iterations, data, config, distance_matrix, print_and_plot=False):
    max_length = max(contig_lengths)
    results = []
    for length in contig_lengths:
        print(f"running for n={length}", end='\r')
        sys.stdout.flush()
        predictions = predict_ages_for_n_contigs(length ,max_length ,iterations, data, config, distance_matrix, how="first")
        results += [predictions]
    if print_and_plot:
        for i in range(len(results)):
            plot_delta_dist(results[i] ,contig_lengths[i])
        plt.legend()
        plt.show()

    return results


def plot_delta_dist(predictions,label):
    predictions.index = [i for i in range(len(predictions[0]))]
#     sns.distplot(predictions[0],label=label)
    print(np.mean(predictions[0]))
    pearson_r = scipy.stats.pearsonr(predictions.iloc[:,0],predictions.iloc[:,1])
#     print([(predictions[0][i],predictions[1][i]) for i in range(len(predictions[0]))])
    print(len(predictions[0]))
    print(f"pearson-r:{pearson_r}")
    delta = [abs(predictions[0][i]-predictions[1][i]) for i in range(len(predictions[0]))]
    print(f"naive mean: {np.mean(delta)}")
    mu, std = scipy.stats.norm.fit(delta)
    print(f"mean erroro:{mu},std:{std}")
    print(f"mean square error:{sklearn.metrics.mean_squared_error(predictions[0], predictions[1])}\n")
    sns.distplot(delta,label=label)

def calc_pearson(predictions):
    return scipy.stats.pearsonr(predictions.iloc[:, 0], predictions.iloc[:, 1])[0]

def calc_spearman(predictions):
    return scipy.stats.spearmanr(predictions.iloc[:, 0], predictions.iloc[:, 1])[0]

def calc_delta_list(predictions):
    return [abs(predictions.iloc[i,0] - predictions.iloc[i,1]) for i in range(len(predictions.iloc[:,0]))]
    # return [abs(predictions.iloc[0,i] - predictions.iloc[1,i]) for i in range(len(predictions.iloc[0]))]

def mean_delta(predictions):
    return np.mean(calc_delta_list(predictions))

def mean_sqaure_error(predictions):
    return sklearn.metrics.mean_squared_error(predictions.iloc[:,0], predictions.iloc[:,1])

def plot_age_prediction(results, rf_results=None, color="b", label="label_placeholder", ax=None, print_stats=False,positions=None):
    if ax is None:
        fig, ax = plt.subplots(2, 2, figsize=(10, 10))
    # if rf_results is not None:
    #     results = [rf_results] + results
    n = len(results)
    if rf_results is None:
        x_values = [i for i in range(1, n+1)]
        x_tick_locations = [i for i in range(n+1)]
        x_ticks = ["RF"] + [str(i) for i in range(1, n+1)]
        print(x_ticks)
        # colors = ["r"] + [color] * (n - 1)
    # else:
        # x_values = [i for i in range(1, n + 1)]
        # x_ticks = [str(i) for i in range(1, n + 1)]
        # colors = [color] * n
    pearson_rs = [calc_pearson(predictions) for predictions in results]
    spearman_rs = [calc_spearman(predictions) for predictions in results]
    mean_errors = [mean_delta(predictions) for predictions in results]
    mean_square_errors = [mean_sqaure_error(predictions) for predictions in results]
    delta_list = [calc_delta_list(predictions) for predictions in results]
    if print_stats:
        for i, deltas in enumerate(delta_list[1:]):
            print(f"1-long vs. {i + 2} long, mann whitney:", scipy.stats.mannwhitneyu(delta_list[0], deltas))
    ax[0, 0].scatter(x_values, pearson_rs, c=color, label=label)
    ax[0, 0].plot(*calc_trendline([i for i in range(1, n)], pearson_rs[1:]), c=color)
    ax[0, 0].set_xticks([i for i in range(n+1)])
    ax[0, 0].set_xticklabels(x_ticks)
    ax[0, 0].legend()
    print("pearson-r between pearson-rs and length:",scipy.stats.pearsonr([i for i in range(1, n + 1)], pearson_rs))
    print(pearson_rs,spearman_rs)

    ax[0, 0].set_title("pearson-r")
    # ax[0, 1].scatter(x_values, mean_errors, c=color)
    # ax[0, 1].plot(*calc_trendline([i for i in range(1, n)], mean_errors[1:]), c=color)
    # ax[0, 1].set_xticks(x_tick_locations)
    # ax[0, 1].set_xticklabels(x_ticks)
    # ax[0, 1].set_title("mean_error")
    ax[0, 1].scatter(x_values, spearman_rs, c=color, label=label)
    ax[0, 1].plot(*calc_trendline([i for i in range(1, n)], spearman_rs[1:]), c=color)
    ax[0, 1].set_xticks([i for i in range(n + 1)])
    ax[0, 1].set_xticklabels(x_ticks)
    ax[0, 1].legend()
    ax[0, 1].set_title("spearman-r")
    ax[1, 1].scatter(x_values, mean_square_errors, c=color)
    ax[1, 1].set_xticks(x_tick_locations)
    ax[1, 1].set_xticklabels(x_ticks)
    ax[1, 1].plot(*calc_trendline([i for i in range(1, n)], mean_square_errors[1:]), c=color)
    ax[1, 1].set_title("mean_square_error")
    if positions is None:
        positions = range(2,n+2)
    ax[1, 0].boxplot(delta_list,positions=positions)
    ax[1, 0].set_xticks(x_tick_locations)
    ax[1, 0].set_xticklabels(x_ticks)

    return ax
#     plt.show()


def plot_age_prediction_compact(results, rf_results=None, color="b", label="label_placeholder", ax=None, print_stats=False,positions=None,title="placeholder"):
    n = len(results)
    x_values = [i for i in range(1, n + 1)]
    x_tick_locations = [i for i in range(1, n + 1)]
    x_ticks = [str(i) for i in range(1, n + 1)]
    # if rf_results is not None:
    #     results = [rf_results] + results
    fig, ax1 = plt.subplots()
    pearson_rs = [calc_pearson(predictions) for predictions in results]
    mean_square_errors = [mean_sqaure_error(predictions) for predictions in results]
    ax1.set_title(title,fontsize=20)
    color = 'tab:blue'
    ax1.set_xlabel('Slice length',fontsize=20)
    ax1.set_ylabel('Pearson-r', color=color,fontsize=20)
    ax1.scatter(x_values, pearson_rs, color=color)
    ax1.plot(x_values, pearson_rs, color=color)
    ax1.set_xticks(x_tick_locations)
    ax1.set_xticklabels(x_ticks)

    ax1.tick_params(axis='y', labelcolor=color)

    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
    color = 'tab:red'
    ax2.set_ylabel('MSE', color=color,fontsize=20)  # we already handled the x-label with ax1
    ax2.scatter(x_values, mean_square_errors, color=color)
    ax2.plot(x_values, mean_square_errors, color=color)
    ax2.tick_params(axis='y', labelcolor=color)

    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    plt.show()



def scatter_predictions(predictions, label):
    plt.scatter(predictions[0], predictions[1], label=label)
    plt.plot(*calc_trendline(predictions[0], predictions[1]))
    plt.show()


def calc_trendline(x, y):
    z = np.polyfit(x, y, 1)
    p = np.poly1d(z)
    return x, p(x)


def predict_ages_with_rf(data_for_rf, samples_to_exclude=[]):
    ages_for_rf = data_for_rf.ages
    ages_for_rf = ages_for_rf.set_index("SampleID")
    ages_for_rf.index = [str(age) for age in ages_for_rf.index.values]
#     ages_for_rf = ages_for_rf.loc[:,["SampleID","Age_at_Collection"]]
    ages_for_rf = ages_for_rf.loc[data_for_rf.ages["SampleID"],"Age_at_Collection"]
    abundance_for_rf = data_for_rf.abundance.T
    merged = abundance_for_rf.join(pd.DataFrame(ages_for_rf))
    merged.drop(samples_to_exclude,axis=0,inplace=True)
#     print(merged["Age_at_Collection"].value_counts())
    regr = RandomForestRegressor(random_state=0, oob_score=True)
    X = merged.iloc[:,:-1]
    y = merged.iloc[:,-1]
    regr.fit(X, y)
    # print("oob_score=:",regr.oob_score_)
    oob_results = pd.DataFrame([regr.oob_prediction_,y])
    return regr, oob_results.T


def add_rf_to_plot(ax,rf_results):
    rf_results = pd.DataFrame(rf_results)
    print(type(rf_results))
    pearson_r = calc_pearson(rf_results)
    mean_error = mean_delta(rf_results)
    mean_square_error = mean_sqaure_error(rf_results)
    delta_list = calc_delta_list(rf_results)
    ax[0,0].scatter([0],[pearson_r],c="r",label="RF")
    ax[0,1].scatter([0],[mean_error],c="r",label="RF")
    ax[1,1].scatter([0],[mean_square_error],c="r",label="RF")
    ax[1,0].boxplot(delta_list)


def rf_ages_for_n_contigs(iterations, data, config, how="first"):
    true = []
    predicted = []
    subjects = list(data.ages["Subject_ID"].unique())
    for i in range(iterations):
        all_predictions = []
        subject1 = random.choice(subjects)
        all_samples = data.sample_ids_by_subject(subject1)
        #         start_index = random.randint(0,len(all_samples)-n)
        start_index = random.randint(0, len(all_samples)-1)
        sample_to_predict = all_samples[start_index]
        age_to_predict = data.get_ages_of_samples([sample_to_predict]).iloc[0]
        print(f"i={i}", end='\r')
        sys.stdout.flush()
        regr = predict_ages_with_rf(data,all_samples)[0]
        prediction = regr.predict(data.abundance[sample_to_predict].values.reshape(1, -1))
        # print(mean_predictions)
        true += [age_to_predict]
        predicted += list(prediction)
    results = pd.DataFrame([true, predicted]).T
    results = results.dropna()
    results.index = [i for i in range(len(results[0]))]
    # print(results)

    return results


def rf_ages_for_n_contigs(iterations, data, config, how="first"):
    true = []
    predicted = []
    subjects = list(data.ages["Subject_ID"].unique())
    for i in range(iterations):
        all_predictions = []
        subject1 = random.choice(subjects)
        all_samples = data.sample_ids_by_subject(subject1)
        #         start_index = random.randint(0,len(all_samples)-n)
        start_index = random.randint(0, len(all_samples)-1)
        sample_to_predict = all_samples[start_index]
        age_to_predict = data.get_ages_of_samples([sample_to_predict]).iloc[0]
        print(f"i={i}", end='\r')
        sys.stdout.flush()
        regr = predict_ages_with_rf(data,all_samples)[0]
        prediction = regr.predict(data.abundance[sample_to_predict].values.reshape(1, -1))
        # print(mean_predictions)
        true += [age_to_predict]
        predicted += list(prediction)
    results = pd.DataFrame([true, predicted]).T
    results = results.dropna()
    results.index = [i for i in range(len(results[0]))]
    # print(results)

    return results

