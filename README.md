# TO DO

* Simulation of activation: might want to create one file where parameters are defined. Instead of having the same code copied into simultation, analysis and plotting.

# Goal of the study
When a researcher wants to combine data from individual fMRI data-analyses, then we consider two options (not exclusive):
 * Image based meta-analysis
 * Third level of a GLM

In this methodological paper, we want to compare these two approaches. We will use **simulated null data** as well as **resting state fMRI data** to calculate point estimates of the true population  (null-) effect and construct 95% confidence intervals (CIs) around these. Then we will calculate:
 * the average standardized bias
 * the average confidence interval length
 * the average coverage of the CIs


### IBMA
This approach is rooted into the meta-analysis literature. For every voxel, we calculate a weighted average using effect size estimates for the true population effect. There are several options for the weights, as well as calculating the CIs. Some of these are discussed in this paper.

### Third level GLM
As is standard in fMRI, statistical modeling in a single study consists of applying a two-stage procedure. First, individual time series are modeled within each subject using a General Linear Model (GLM) in each voxel. Then several subjects are combined in the second stage using again a GLM. The idea is to pull this through at a third stage in which studies are combined using a GLM. Intersubject, between-subject and between-study variability can be accounted for in each stage.

# Approach
The idea is to use fMRI null data to calculate the properties mentioned above. Indeed using:

* Simulations in R
* Resting state fMRI data from the 1000 functional connectome project (as in Eklund et al., 2016)
* Resting state fMRI data from the 1000 functional connectomes project (as in Eklund et al., 2016)
* Possibly an approach provided by Slotnick (2017) in which task based fMRI data is taken. But odd trials are contrasted with even trials from the same task. This creates a null data situation as well.

The reason for this extended approach is to provide as much validation as possible. Each method has its advantages and disadvantages.


# Structure of repository
* [1_Scripts](https://github.com/NeuroStat/SimulationGit/tree/master/1_Scripts) contains the scripts.
* [2_Analyses](https://github.com/NeuroStat/SimulationGit/tree/master/2_Analyses) contains the R code for analyses.
* [3_Reports](https://github.com/NeuroStat/SimulationGit/tree/master/3_Reports) contains the latest report only.


# Link with OSF
This project is hosted on the [Open Science Framework][1].
 [1]: https://osf.io/t92bd/

# Getting started with Git
* https://git-scm.com/doc
* https://www.youtube.com/user/GitHubGuides/videos

# Contact
[Webpage][Han Bossier] or [e-mail](mailto: Han.Bossier@Ugent.be).

[Han Bossier]: http://telefoonboek.ugent.be/nl/people/802001626303
[Freya Acar]: https://telefoonboek.ugent.be/nl/people/802001860820
[Ruth Seurinck]: http://telefoonboek.ugent.be/nl/people/801001629152
[Beatrijs Moerkerke]: http://telefoonboek.ugent.be/nl/people/801001453542

# References
Eklund A, Nichols TE, Knutsson H, Cluster failure: Why fMRI inferences for spatial extent have inflated false-positive rates., (2016), "Proc Natl Acad Sci USA. doi: [10.1073/pnas.1602413113](http://www.pnas.org/content/113/28/7900.abstract?sid=1f6fd91a-988c-4a69-80e4-a19ae1951614)


Slotnick S.D., Resting-state fMRI data reflects default network activity rather than null data: A defense of commonly employed methods to correct for multiple comparisons., 2017, Cogn Neurosci, doi: [10.1080/17588928.2016.1273892](http://www.tandfonline.com/doi/full/10.1080/17588928.2016.1273892?scroll=top&needAccess=true).



First Created: 12/01/16
