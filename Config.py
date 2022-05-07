
class Config:

    def __init__(self,abundance_threshold = 0, samples_above_abundance_threshold=2, max_age = 3000,min_max_age=0,window_size = 30,top_predictors_for_age_analysis=0,k=30):
        self.abundance_threshold = abundance_threshold
        self.max_age = max_age #maximal age to include in dataset
        self.min_max_age=min_max_age #keep only subject who's maximal age greater than this parameter
        self.samples_above_abundance_threshold = samples_above_abundance_threshold
        # Interpolation Parameters - fill one and leave the other as None
        self.number_of_days = None
        self.day_interval = 10
        self.window_size = window_size
        self.k_day_bins = True
        self.k = k


        #change name of "time-threshold" to something that makes more sense (scaling factor? norm factor?)

        #Distance Calculation Parameters
        # normalization methods - "regular" - substract minimum and divide by maximum.
        # "log" - take log and clip to [0,1]. "linear" - divide by threshold and clip
        self.metric = "braycurtis" #"braycurtis" or "weighted_unifrac"
        self.time_prior_normalization = "linear"
        self.time_threshold = 300
        self.log_base = None
        self.time_prior_weight = 0 #number in [0,1]

        #change name of local alignment?
        #add options for choosing best local alignment based on maximum length (with distance threshold) or best score (with length minimum)

        # Alignment parameters
        self.calc_local_alignment = True
        self.calc_global_alignment = True
        self.global_step_pattern="symmetric2" # hard coded for now
        self.align_threshold_value = None #for local alignment
        self.align_threshold_quantile = 0.1 #for local alignment

        #Output settings:
        self.number_of_top_hitters = 5
        self.number_of_bottom_hitters = 5
        self.show_pairwise_pcoa = True
        self.show_tree_pcoa = False
        self.show_interpolation_pcoa = True

        # Age analysis parameters
        self.top_predictors_for_age_analysis = top_predictors_for_age_analysis