# Between-study variability

# Important Note

Data and analyses are incorporated in the overview paper on effect sizes. Hence, sections of this folder will be copied and modified to
https://github.com/NeuroStat/ESfMRI

# Introduction

Section to try and estimate between-study variability in neuroimaging (focussed on fMRI).

Data is obtained by searching for *"pain"* (at 11/01/2017) in the [neurovault](www.neurovault.org) database and manually checking all results.

# Topic
The idea is to search and gather full SPM's on *Pain* versus *No pain*.
We want to obtain the distribution of values for between-study heterogeneity in areas involved with experiencing pain. To get these areas, we create a mask using [neurosynth](www.neurosynth.org) searching for the term **pain**.
We obtain an automated meta-analysis of 420 studies. From here we use reverse inference with a FDR at 0.01, the standard neurosynth procedure to create a mask. Reverse inference is equal to: P(Term|Activation)<sup>1</sup>.
This is the map called *pain_pFgA_z_FDR_0.01.nii*. After estimating the between-study heterogeneity, we will use this mask to obtain the distribution.

# Studies

The data is stored in the *Data/* folder and the associated papers in *Papers/*. However, these folders are not pushed to the Github repository as they are too large. <br>
The database contains the following studies:


| Study        | Sample size           | Type  |  Contrast |
| ------------- |:-------------|:-----:|:-----:|
|Braboszcz 2017 | 17            | T-map | Painful images >Painless images (Normal state) |
| Hebestreit 2017      | 23      |   T-map | Ammonia stimulation (trigeminal pain) > Air in both sessions (medication and placebo) |
| Tamm 2017 | 86      |  T-map | Pain > no pain [no covariates] |
| Karjalainen 2017 | 35      |  T-map | Main effect of vicarious pain |
| Atlas 2010 | 15      |  beta-map<sup>2</sup> | Thermal high vs low stimulation pain |
| Wager 2013 | 15      |  beta-map<sup>2</sup> | Somatic pain vs baseline |
| Kano 2017 | 15      |  beta-map<sup>2</sup> | Visceral pain vs baseline |
| Rubio 2015 | 15      |  beta-map<sup>2</sup> | Visceral pain vs baseline |
| Unpublished | 15      |  beta-map<sup>2</sup> | Mechanical high pressure pain vs baseline |
| Unpublished | 15      |  beta-map<sup>2</sup> | Mechanical medium pressure pain vs baseline |
| Patil 2017 | 46      |  T-map | Outcome of pain induced to others versus baseline |
| Maumet 2016 | Total = 334<sup>3</sup>      |  T-maps | Pain versus baseline |


* Note that we only include studies doing whole-brain analyses.
* <sup>2</sup>Data comes from Kragel et al, 2017. We need to pool the subjects using OLS. I assume the provided data from this study are first level beta parameter estimates each subject. Should verify with Tor Wager.
* <sup>3</sup>Study of Maumet et al., 2016 is about the *nidm* data structure. However, it contains **21** studies about pain. Note that these come from the same site! Might need to control for when estimating between-study variability.

## Notes
In the following study: https://neurovault.org/collections/504/, it seems we can create beta maps for high vs low pain. Using OLS, it might be possible to have a study with approx. 30 subjects.

# Estimation


<sup>1</sup> http://neurosynth.org/faq/#q15
