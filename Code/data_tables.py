from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import Config

class ta_data():

    def __init__(self,metadata,config,counts=None,tree=None,abundance=None):
        self.config = config
        if not isinstance(abundance,type(None)):
            self.abundance = abundance #samples in columns, taxa in row. RELATIVE ABUNDANCE ONLY
        else:
            self.abundance = counts.div(counts.sum(axis=0), axis=1)
        if config.abundance_threshold>0:
            self.abundance = self.abundance[(self.abundance > config.abundance_threshold).sum(axis=1) >= config.samples_above_abundance_threshold]
        self.metadata = metadata
        self.filter_min_max_age()
        samples_to_keep = np.intersect1d(self.metadata["SampleID"], self.abundance.columns)
        # print(f"abundance before:{self.abundance.shape}, metadata before:{metadata.shape}")
        self.metadata = self.metadata[self.metadata["SampleID"].isin(samples_to_keep)]
        self.abundance = self.abundance[samples_to_keep]
        # print(f"abundance after:{self.abundance.shape}, metadata after:{self.metadata.shape}")
        if not all([abs(sum(self.abundance.loc[:,column])-1)<0.15 for column in self.abundance.columns]):
            print("warning - abundance dowsn't sum up to 1")
        #making sure that the sum of each column is 1
        self.counts=counts
        self.tree=tree
        assert "Age_at_Collection" in self.metadata.columns
        assert "SampleID" in self.metadata.columns
        assert "Subject_ID" in self.metadata.columns
        self.ages = self.metadata.sort_values(by=['Subject_ID', "Age_at_Collection"])
        # samples_to_drop = [sample for sample in self.ages.loc[:,"SampleID"] if sample not in self.abundance.columns]
        # print("ages before drop",self.ages.shape,"abundance before drop",self.abundance.shape)
        # self.ages=self.ages[~self.ages["SampleID"].isin(samples_to_drop)]
        self.ages=self.ages.loc[:,["Subject_ID","SampleID","Age_at_Collection"]]
        self.phyla = None

    def copy(self):
        new = ta_data(metadata=self.metadata.copy(),abundance=self.abundance.copy(),config=self.config)
        new.ages = self.ages
        return new

    def filter_min_max_age(self):

        subject_filter = self.metadata.groupby("Subject_ID").max()["Age_at_Collection"] > self.config.min_max_age
        filtered_metadata = self.metadata[self.metadata["Subject_ID"].isin(list(subject_filter[subject_filter == True].index))]
        filtered_metadata = filtered_metadata[filtered_metadata["Age_at_Collection"] < self.config.max_age]
        self.metadata = filtered_metadata


    def sample_ids_by_subject(self, subject):
        # returns all the sample ID's of the subject, excluding those missing from the OTU table
        ids = self.ages[self.ages["Subject_ID"] == subject]["SampleID"]
        ids = [id for id in ids if id in self.abundance.columns]

        return ids

    def subject_of_sample(self,sample):
        return self.ages[self.ages["SampleID"] == sample]["Subject_ID"].iloc[0]

    def abundance_by_subject(self,subject):
        ids = self.sample_ids_by_subject(subject)
        return self.abundance[ids]

    def all_subjects(self):
        return list(self.ages["Subject_ID"].unique())


    def counts_by_subject(self,subject):
        ids = self.sample_ids_by_subject(subject)
        return self.counts[ids]

    def ages_by_subject(self,subject):
        metadata=self.ages
        ages_by_subject = metadata[metadata["Subject_ID"] == subject]["Age_at_Collection"]
        ages_by_subject.index = metadata[metadata["Subject_ID"] == subject]["SampleID"]
        return ages_by_subject

    def age_of_sample(self, sample):
        return self.ages[self.ages["SampleID"] == sample]["Age_at_Collection"]

    def get_ages_of_samples(self,samples):
        # ages = self.ages[self.ages["SampleID"].isin(samples)]["Age_at_Collection"]
        # ages.index = samples
        ages =  [int(self.age_of_sample(sample)) for sample in samples]
        return pd.DataFrame(ages,index=samples).iloc[:,0]

    def get_subjects_of_samples(self,samples):
        # ages = self.ages[self.ages["SampleID"].isin(samples)]["Age_at_Collection"]
        # ages.index = samples

        return [self.subject_of_sample(sample) for sample in samples]


    def add_samples_and_ages(self,subject,samples,ages):
        pass

    def get_sample(self,sample):
        return self.abundance[sample]



    def plot_age_ranges(self):
        plt.style.use("ggplot")
        subjects = list(self.ages["Subject_ID"].unique())
        x = []
        y = []
        for i, subject in enumerate(subjects):
            ages = list(self.ages_by_subject(subject))
            #     ages = [age for age in ages if age<800]
            n = len(ages)
            y += [i] * n
            x += ages

        plt.rcParams["figure.figsize"] = (15, 15)
        plt.xlabel("Age (days)",fontsize=20)
        plt.ylabel("Subject-ID", fontsize=20)
        plt.yticks(range(len(subjects)), subjects)
        plt.scatter(x, y)
        # plt.yticks(rotation=10)
        for i, subject in enumerate(subjects):
            ages = list(self.ages_by_subject(subject))
            #     ages = [age for age in ages if age<800]
            x_coor = [min(ages), max(ages)]
            y_coor = [i, i]
            plt.plot(x_coor, y_coor)

        plt.show()

def merge(data_table_1, data_table_2):
    abundance_all = pd.concat([data_table_1.abundance, data_table_2.abundance], axis=1, join="inner")
    ages = pd.concat([data_table_1.ages,data_table_2.ages],axis=0)
    config = Config.Config(abundance_threshold=0)
    return ta_data(metadata=ages,abundance=abundance_all,config=config)
