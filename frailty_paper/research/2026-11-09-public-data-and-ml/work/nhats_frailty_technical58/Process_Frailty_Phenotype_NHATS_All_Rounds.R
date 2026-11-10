# R Code for deriving Frailty Phenotype from Multiple NHATS Rounds
#
# Prior to running this code:
# (1) Please take note of your working directory. The paths specified below
# should be relative to this working directory.
#
# (2) Please make sure appropriate data has been downloaded. See section 5.1
# "Downloading the Data".
#
# (3) Be sure to install the following required packages in R:
#	tidyverse (for general data processing)
#	haven (for reading STATA or SAS files)

################################################################################
# User Input Section

# Relative path from your current R working directory to the parent
# directory which contains data files for all rounds.
data_parent_directory <- "Data_Alternative"

# Select whether cutpoints should be pulled from the Bandeen-Roche (2015) paper
# or Fried (2001)
cutpoint_type = "NHATS"
#cutpoint_type = "CHS"

# Define the relative path from your current R working directory to the name of
# the directory in which you would like the derived datasets to be deposited.
data_output_directory <- "Data_Derived"

# Set type of file output: "rdata" or "csv". The "rdata" format will retain
# variable and value labels.
#data_output_format <- "csv"
data_output_format <- "rdata"

################################################################################
# Data Processing Section (edit at your own risk!)

library(tidyverse)

################################################################################
# (1) Cutpoint Selection
################################################################################
if(cutpoint_type == "NHATS"){
    # Low BMI Cutpoint
    bmi_lower_bound = 18.5

    # Weak Grip Parameters
    # sex-specific bmi cutoffs for determining strength cutoff
    f_bmi = c(18.5, 25, 30)
    m_bmi = c(18.5, 25, 30)
    # sex/bmi-specific grip cutoffs
    f_grip = c(12.66, 12.66, 12.92, 12.23)
    m_grip = c(22.75, 22.75, 25.59, 25.96)

    # Slow Gait Parameters
    # sex-specific height cutoff for using higher/lower gait speed cutoff
    f_height = 1.59
    m_height = 1.73

    # sex/height-specific gait speed cutoffs
    f_gait = c(0.35, 0.45)
    m_gait = c(0.5, 0.58)

}else if(cutpoint_type == "CHS"){
	# Low BMI Cutpoint
	bmi_lower_bound = 18.5

	# Weak Grip Parameters
	f_bmi = c(23, 26, 29)
	m_bmi = c(24, 26, 28)
	f_grip = c(17, 17.3, 18, 21)
	m_grip = c(29, 30, 30, 32)

	# Slow Gait Parameters

	# sex-specific height cutoff for using higher/lower gait speed cutoff
	f_height = 1.59
	m_height = 1.73

	# sex-specific gait speed cutoff for taller and shorter group
	f_gait = c(0.65, 0.76)
	m_gait = c(0.65, 0.76)
}

################################################################################
# (2) File Identification, Pre-Processing
################################################################################
df <- tibble(
	file = dir(
		data_parent_directory,
		full.names = T,
		recursive = T,
		ignore.case = T,
		pattern = "NHATS_Round_[0-9].*_SP_File.*\\.dta"
	)
) %>% mutate(
	round = as.numeric(gsub(".*Round_","",gsub("_SP_.*","",file)))
) %>% filter(
	(!is.na(round))
) %>% mutate(
	data = map2(
		file, round,
		~haven::read_dta(
			.x,
			col_select = c(
				starts_with("spid"),
				ends_with(paste0("r",.y, "dresid")),

				# demographics
				ends_with("intvrage"),
				ends_with("dgender"),

				# anthropometry components
				ends_with("howtallft"),
				ends_with("howtallin"),
				ends_with("currweigh"),
				ends_with("lst10pnds"),
				ends_with("trytolose"),

				# low energy items
				ends_with("loenlmtat"),
				ends_with("lowenergy"),

				# low activity items
				ends_with("vigoractv"),
				ends_with("evrgowalk"),

				# mobility items: determines exclusion for gait assessment
				ends_with("insdwlchr"),
				ends_with("insdsctr"),

				# performance activity elligibility items: determines exclusion for gait assessment
				ends_with("wlkdsself"),

				# gait speed items
				ends_with("wlkc1secs"),
				ends_with("wlk1hndr"),
				ends_with("wlkc2secs"),
				ends_with("wlk2hndr"),
				ends_with("rsn31wkc"),
				ends_with("rsn21wkc"),
				ends_with("rsn11wkc"),
				ends_with("dwlkadm"),
				ends_with("wlkc1rslt"),

				# grip strength items
				ends_with("surghdwrt"),
				ends_with("surgyside"),
				ends_with("surgarmsh"),
				ends_with("sidsurgar"),
				ends_with("flruppain"),
				ends_with("sideflrup"),
				ends_with("grp1rdng"),
				ends_with("grp2rdng"),
				ends_with("dgripadm"),
				ends_with("grp1reslt"),
				ends_with("grp1noat1"),
				ends_with("grp1noat2"),
				ends_with("grp1noat3")
			)
		) %>% haven::zap_labels() %>% rename_with(
		    .fn = function(v) sub("r[0-9]+dgender","r__dgender", v)
		) %>% rename_with(
			.fn = function(v) sub(.y,"__", v)
		)
	)
) %>% unnest(
	data
) %>% group_by(
	spid
################################################################################
# (3) Derive intermediate Variables
################################################################################
) %>% mutate(
    # sex was not reported the same way during round 2 so we need to
    # aggregate across rounds
	female = ifelse(
	    sum(r__dgender==1,na.rm=T)>sum(r__dgender==2,na.rm=T),
	    0,
	    1
	),
	.after = "spid"
) %>% ungroup %>% mutate(
    frailty_inapplicable = 1*(
        (r__dresid %in% c(3,5,6,7,8))|
            ((round==1)&(r__dresid == 4))
    ),
    # anthropometrics
    height_m = ifelse(
        (hw__howtallft>0)&(hw__howtallin>=0),
        (hw__howtallft*12+hw__howtallin)/39.37,
        NA
    ),
    weight_kg = ifelse(
        hw__currweigh > 0,
        hw__currweigh/2.205,
        NA
    ),
    bmi = weight_kg/height_m^2,

	# Weight Reduction and Low BMI
	low_bmi = 1*(bmi <= bmi_lower_bound),
	weight_lower = ifelse(
		(hw__lst10pnds < 1),
		NA,
		ifelse(
			(hw__trytolose == 1) | (hw__lst10pnds == 2),
			0,
			# were trying, or did not lose
			ifelse(
				hw__trytolose == 2,
				# was not trying
				1,
				NA
			)
		)
	),

	# speed
	gait1 = ifelse(
		(wa__wlkc1secs > 0)&(wa__wlk1hndr > -1),
		3/(wa__wlkc1secs + wa__wlk1hndr/100),
		NA

	),
	gait2 = ifelse(
		(wa__wlkc2secs > 0)&(wa__wlk2hndr > -1),
		3/(wa__wlkc2secs+wa__wlk2hndr/100),
		NA
	),
	gait_speed_max = pmax(gait1,gait2,na.rm = T),
	gait_unable_unsafe = ifelse(
	    # we will only use this indicator if there is no recorded gait speed
		is.na(gait_speed_max)&(
		    # participant always uses wheelchair OR scooter OR
		    (mo__insdsctr==1)|(mo__insdwlchr==1)|
		    # indicates cannot walk short distance by self OR
		    (pe__wlkdsself==2)|
		    # someone indicates safety concern OR
		    (wa__rsn31wkc == 1)|(wa__rsn21wkc == 1)|(wa__rsn11wkc == 1)|
		    # they do not complete the first trial
		    (wa__wlkc1rslt == 2)
		), 1, NA
	),

	# grip strength
	grip_strength1 = ifelse(
		gr__grp1rdng > 0, gr__grp1rdng, NA
	),
	grip_strength2 = ifelse(
		gr__grp2rdng > 0, gr__grp2rdng, NA
	),
	grip_strength_max = pmax(grip_strength1,grip_strength2, na.rm = T),
	grip_unable_unsafe = ifelse(
		is.na(grip_strength_max)&(
		    (gr__grp1reslt == 2)|
		    # first trial attempted but not completed
			(gr__dgripadm == 4)|
			# Not Eligible due to exclusion criteria
			(gr__grp1noat1 == 1)|(gr__grp1noat2 == 1)|(gr__grp1noat3 == 1)
			# at least one person felt unsafe
		), 1, NA
	),
################################################################################
# (4) Frailty Component Derivation
################################################################################
    # low energy
    low_energy_component = ifelse(
        (ss__loenlmtat == 1)&(ss__lowenergy==1),
        1,
        ifelse(
            (ss__lowenergy == 2) | (ss__loenlmtat == 2),
            0, NA
        )
    ),

    # low activity
    low_activity_component = ifelse(
        (pa__vigoractv == 2)&(pa__evrgowalk == 2),
        1,
        ifelse(
            (pa__vigoractv < 0)&(pa__evrgowalk < 0),
            NA, 0
        )
    ),

    # weight loss (combined)
    weight_loss_component = 1*((low_bmi == 1)|(weight_lower == 1)),

    # slow gait
    slow_gait_component = NA,
    # use strict cutpoints when height and gait are measured
    slow_gait_component = ifelse(
        (!is.na(height_m))&(!is.na(gait_speed_max))&(!is.na(female)),
        ((female == 1)&(height_m<=f_height)&(gait_speed_max <= f_gait[1]))|
            ((female == 1)&(height_m>f_height)&(gait_speed_max <=f_gait[2]))|
            ((female == 0)&(height_m<=m_height)&(gait_speed_max <=m_gait[1]))|
            ((female == 0)&(height_m>m_height)&(gait_speed_max <=m_gait[2])),
        slow_gait_component
    ),
    # use max/min cutpoints when height is not measured
    slow_gait_component = ifelse(
        (is.na(slow_gait_component))&
        (is.na(height_m))&
        (!is.na(gait_speed_max))&
        (!is.na(female))&(
            ((female == 1)&(gait_speed_max <= min(f_gait)))|
                ((female == 0)&(gait_speed_max <= min(m_gait)))
        ), 1, slow_gait_component
    ),
    slow_gait_component = ifelse(
        (is.na(slow_gait_component))&
        (is.na(height_m))&
        (!is.na(gait_speed_max))&
        (!is.na(female))&(
            ((female == 1)&(gait_speed_max >= max(f_gait)))|
            ((female == 0)&(gait_speed_max >= max(m_gait)))
        ), 0, slow_gait_component
    ),
    slow_gait_component = ifelse(
        (is.na(slow_gait_component))&
        (!is.na(gait_unable_unsafe))&
        (gait_unable_unsafe == 1),
        1,
        slow_gait_component
    ),

    # weak_ grip
    # use strict cutpoints when BMI and grip are measured
	weak_grip_component = ifelse(
		(!is.na(female))&(!is.na(bmi))&(!is.na(grip_strength_max)),
		((female == 1)&(bmi<=f_bmi[1])&(grip_strength_max <= f_grip[1]))|
			((female == 1)&(bmi>f_bmi[1])&(bmi<=f_bmi[2])&(grip_strength_max <=f_grip[2]))|
			((female == 1)&(bmi>f_bmi[2])&(bmi<=f_bmi[3])&(grip_strength_max <=f_grip[3]))|
			((female == 1)&(bmi>f_bmi[3])&(grip_strength_max <=f_grip[4]))|
			((female == 0)&(bmi<=m_bmi[1])&(grip_strength_max <=m_grip[1]))|
			((female == 0)&(bmi>m_bmi[1])&(bmi<=m_bmi[2])&(grip_strength_max <=m_grip[2]))|
			((female == 0)&(bmi>m_bmi[2])&(bmi<=m_bmi[3])&(grip_strength_max <=m_grip[3]))|
			((female == 0)&(bmi>m_bmi[3])&(grip_strength_max <= m_grip[4])),
		NA
	),
    # use max/min cutpoints when BMI not measured and grip is measured
	weak_grip_component = ifelse(
		(is.na(weak_grip_component)) &
		(!is.na(female)) &
		(is.na(bmi)) &
		(!is.na(grip_strength_max))&(
			((female == 1)&(grip_strength_max <=min(f_grip)))|
			((female == 0)&(grip_strength_max <=min(m_grip)))
		), 1, weak_grip_component
	),
    weak_grip_component = ifelse(
        is.na(weak_grip_component) &
        (!is.na(female)) &
        (is.na(bmi)) &
        (!is.na(grip_strength_max))&(
            ((female == 1)&(grip_strength_max >max(f_grip)))|
            ((female == 0)&(grip_strength_max >max(m_grip)))
        ), 0, weak_grip_component
    ),

    # use unable/unsafe criteria otherwise
	weak_grip_component = ifelse(
	    is.na(weak_grip_component) &
	        (!is.na(grip_unable_unsafe))&
	        (grip_unable_unsafe == 1),
		1,
		weak_grip_component
	),

    # Weak grip reason
    weak_grip_reason = ifelse(
        frailty_inapplicable == 1,
        -1,
        NA
    ),
    weak_grip_reason = ifelse(
        is.na(weak_grip_component)&is.na(weak_grip_reason),
        -9,
        weak_grip_reason
    ),
    weak_grip_reason = ifelse(
        (weak_grip_component==0)&is.na(weak_grip_reason)&!is.na(weak_grip_component),
        -8,
        weak_grip_reason
    ),
    weak_grip_reason = ifelse(
        is.na(weak_grip_reason)&(!is.na(grip_strength_max)),
        6,
        weak_grip_reason
    ),
    weak_grip_reason = ifelse(
        is.na(weak_grip_reason)&(gr__dgripadm == 4)&((pe__surghdwrt ==1)|(pe__surgarmsh ==1)),
        1,
        weak_grip_reason
    ),
    weak_grip_reason = ifelse(
        is.na(weak_grip_reason)&(gr__dgripadm == 4)&(pe__flruppain ==1),
        2, weak_grip_reason
    ),
    weak_grip_reason = ifelse(
        is.na(weak_grip_reason)&(gr__dgripadm == 4),
        5, weak_grip_reason
    ),
    weak_grip_reason = ifelse(
        is.na(weak_grip_reason) & ((gr__grp1noat1 == 1)|(gr__grp1noat2 == 1)|(gr__grp1noat3 == 1)),
        3,
        weak_grip_reason
    ),
    weak_grip_reason = ifelse(
        is.na(weak_grip_reason) & (gr__grp1reslt == 2),
        4,
        weak_grip_reason
    ),

    # slow gait reason
    slow_gait_reason = ifelse(
        frailty_inapplicable == 1,
        -1,
        NA
    ),
    slow_gait_reason = ifelse(
        is.na(slow_gait_component)&is.na(slow_gait_reason),
        -9,
        slow_gait_reason
    ),
    slow_gait_reason = ifelse(
        is.na(slow_gait_reason)&(slow_gait_component==0)&(!is.na(slow_gait_component)),
        -8,
        slow_gait_reason
    ),

    slow_gait_reason = ifelse(
        is.na(slow_gait_reason)&(!is.na(gait_speed_max)),
        5,
        slow_gait_reason
    ),
    slow_gait_reason = ifelse(
        is.na(slow_gait_reason)&((mo__insdsctr==1)|(mo__insdwlchr==1)),
        1,
        slow_gait_reason
    ),
    slow_gait_reason = ifelse(
        is.na(slow_gait_reason) & (pe__wlkdsself == 2),
        2,
        slow_gait_reason
    ),
    slow_gait_reason = ifelse(
        is.na(slow_gait_reason) & ((wa__rsn31wkc == 1)|(wa__rsn21wkc == 1)|(wa__rsn11wkc == 1)),
        3,
        slow_gait_reason
    ),
    slow_gait_reason = ifelse(
        is.na(slow_gait_reason) & (wa__wlkc1rslt == 2),
        4,
        slow_gait_reason
    ),

################################################################################
# (5) Frailty Phenotype Derivation
################################################################################
	# Count Frailty Components
	# the minimum possible score is derived by assuming NAs are all 0
	min_score =
		ifelse(is.na(weak_grip_component),0,weak_grip_component) +
		ifelse(is.na(slow_gait_component),0,slow_gait_component) +
		ifelse(is.na(low_energy_component),0,low_energy_component) +
		ifelse(is.na(low_activity_component),0,low_activity_component) +
		ifelse(is.na(weight_loss_component),0,weight_loss_component),

	# the maximum possible score is derived by assuming NAs are all 1
	max_score =
		ifelse(is.na(weak_grip_component),1,weak_grip_component) +
		ifelse(is.na(slow_gait_component),1,slow_gait_component) +
		ifelse(is.na(low_energy_component),1,low_energy_component) +
		ifelse(is.na(low_activity_component),1,low_activity_component) +
		ifelse(is.na(weight_loss_component),1,weight_loss_component),

	# a frailty phenotype is derived when the range of possible scores falls
	# definitively within the criteria for a given category

    frailty_phenotype = ifelse(
        frailty_inapplicable == 1,
        -1,
        NA
    ),
    frailty_phenotype = ifelse(
        (min_score >=3)&is.na(frailty_phenotype),
        1,
        frailty_phenotype
    ),
	frailty_phenotype = ifelse(
		(min_score > 0)&(max_score<3)&is.na(frailty_phenotype),
		2,
		frailty_phenotype
	),
    frailty_phenotype = ifelse(
        (max_score == 0)&is.na(frailty_phenotype),
        3,
        frailty_phenotype
    ),
    frailty_phenotype = ifelse(
        (min_score>0)&is.na(frailty_phenotype),
        4,
        frailty_phenotype
    ),
    frailty_phenotype = ifelse(
        (max_score<3)&is.na(frailty_phenotype),
        5,
        frailty_phenotype
    ),
    frailty_phenotype = ifelse(
        (min_score==0)&(max_score==5)&is.na(frailty_phenotype),
        -9,
        frailty_phenotype
    ),
    frailty_phenotype = ifelse(
        is.na(frailty_phenotype),
        -8,
        frailty_phenotype
    )
)

################################################################################
# (6) Formatting
################################################################################

df <- df %>% relocate(
    female, .after = "frailty_inapplicable"
) %>% mutate(
    across(
        female:slow_gait_reason,
        .fns = function(x){
            ifelse(
                (frailty_inapplicable==1)&is.na(x),
                -1,
                ifelse(
                    is.na(x),
                    -9,
                    x
                )
            )
        }
    ),
    min_score = ifelse(frailty_inapplicable, -1, min_score),
    max_score = ifelse(frailty_inapplicable, -1, max_score)
) %>% rename_with(
    .cols = frailty_inapplicable:frailty_phenotype,
    .fn = function(v) paste0("r__d_",v)
)

# Assign label to variables
attr(df$r__d_frailty_inapplicable,'label') <-
    "Frailty Inapplicable: No SP interview."
attr(df$r__d_female,'label') <-
    "Female (vs. Male)"
attr(df$r__d_height_m,'label') <-
    "Height (m)"
attr(df$r__d_weight_kg,'label') <-
    "Weight (kg)"
attr(df$r__d_bmi,'label') <-
    "Body Mass Index (kg/m^2)"
attr(df$r__d_low_bmi,'label') <-
    "SP had low derived BMI."
attr(df$r__d_weight_lower,'label') <-
    "SP reported unintentional weight loss."
attr(df$r__d_gait1,'label') <-
    "First Gait Speed Trial"
attr(df$r__d_gait2,'label') <-
    "Second Gait Speed Trial"
attr(df$r__d_gait_speed_max,'label') <-
    "Maximum of Observed Gait Speed Trials"
attr(df$r__d_gait_unable_unsafe,'label') <-
    "Gait speed could not be derived but the SP met some other criteria for slow gait."
attr(df$r__d_grip_strength1,'label') <-
    "First Grip Strength Trial"
attr(df$r__d_grip_strength2,'label') <-
    "Second Grip Strength Trial"
attr(df$r__d_grip_strength_max,'label') <-
    "Maximum of Observed Grip Strength Trials"
attr(df$r__d_grip_unable_unsafe,'label') <-
    "Grip strength could not be derived but the SP met some other criteria for weak grip."
attr(df$r__d_low_energy_component,'label') <-
    "Low Energy component of Frailty Phenotype."
attr(df$r__d_low_activity_component,'label') <-
    "Low Activity component of Frailty Phenotype."
attr(df$r__d_weight_loss_component,'label') <-
    "Weight Loss component of Frailty Phenotype."
attr(df$r__d_slow_gait_component,'label') <-
    "Slow Gait component of Frailty Phenotype."
attr(df$r__d_weak_grip_component,'label') <-
    "Weak Grip component of Frailty Phenotype."
attr(df$r__d_min_score,'label') <-
    "The number of positive and observed frailty criteria, not counting missing criteria. (0 to 5)"
attr(df$r__d_max_score,'label') <-
    "The number of positive or missing frailty criteria. (0 to 5)"
attr(df$r__d_frailty_phenotype,'label') <-
    "Frailty Phenotype or indefinite classification"
attr(df$r__d_weak_grip_reason,'label') <-
    "Trigger for weak grip classification."
attr(df$r__d_slow_gait_reason,'label') <-
    "Trigger for slow gait classification."

# Assign labels to values
yes_no_label <- c(
    "Inapplicable" = -1,
    "No" = 0,
    "Yes" = 1
)

attr(df$r__d_frailty_inapplicable,'labels') <- yes_no_label
attr(df$r__d_female,'labels') <- yes_no_label
attr(df$r__d_low_bmi,'labels') <- yes_no_label
attr(df$r__d_weight_lower,'labels') <- yes_no_label
attr(df$r__d_gait_unable_unsafe,'labels') <- c(
    "Inapplicable" = -1,
    "Yes" = 1
)
attr(df$r__d_grip_unable_unsafe,'labels') <- c(
    "Inapplicable" = -1,
    "Yes" = 1
)
attr(df$r__d_low_energy_component,'labels') <- yes_no_label
attr(df$r__d_low_activity_component,'labels') <- yes_no_label
attr(df$r__d_weight_loss_component,'labels') <- yes_no_label
attr(df$r__d_slow_gait_component,'labels') <- yes_no_label
attr(df$r__d_weak_grip_component,'labels') <- yes_no_label
attr(df$r__d_frailty_phenotype,'labels') <- c(
    "1 Frail" = 1,
    "2 Prefrail" = 2,
    "3 Robust" = 3,
    "4 Not Robust" = 4,
    "5 Not Frail" = 5,
    "-1 Inapplicable" = -1,
    "-8 Indeterminate" = -8,
    "-9 Not Assessed" = -9
)
attr(df$r__d_weak_grip_reason,'labels') <- c(
    "-9 Missing Weakness Component" = -9,
    "-8 Not Weak" = -8,
    "-1 Inapplicable" = -1,
    "1 Weak, Surgery or Injury Exclusion" = 1,
    "2 Weak, Pain Exclusion" = 2,
    "3 Weak, Safety Objection" = 3,
    "4 Weak, Did Not Complete" = 4,
    "5 Weak, Other Reason" = 5,
    "6 Weak, Grip Strength Below Threshold" = 6
)
attr(df$r__d_slow_gait_reason,'labels') <- c(
    "-9 Missing Slowness Component" = -9,
    "-8 Not Slow" = -8,
    "-1 Inapplicable" = -1,
    "1 Slow, Scooter or Wheelchair Use" = 1,
    "2 Slow, Unable to Walk by Self" = 2,
    "3 Slow, Safety Objection" = 3,
    "4 Slow, Did Not Complete" = 4,
    "5 Slow, Gait Speed Below Threshold" = 5
)

################################################################################
### Export CSV File

# create output directory if needed
if(!dir.exists(data_output_directory)){
	dir.create(data_output_directory, recursive = T)
}

# create timestamp to append to data
data_timestamp = format(now(), format = "%Y%m%d_%H%M%S")

if(data_output_format == "csv"){
    write_csv(
        df,
        file.path(
            data_output_directory,
            paste0(
                "NHATS_Frailty_Longitudinal_",
                data_timestamp,
                ".csv"
            )
        )
    )
} else{
    saveRDS(
        df,
        file.path(
            data_output_directory,
            paste0(
                "NHATS_Frailty_Longitudinal_",
                data_timestamp,
                ".rdata"
            )
        )
    )
}


