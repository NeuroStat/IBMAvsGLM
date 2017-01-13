# Goal of the study
When a researcher wants to combine data from individual fMRI data-analyses, then we consider two options (not exclusive):
 * Image based meta-analysis
 * Third level of a GLM

In this methodological paper, we want to compare these two approaches. We will use **simulated null data** as well as **resting state fMRI data** to calculate point estimates of the true population  (null-) effect and construct 95% confidence intervals (CIs) around these. Then we will calculate:
 * the average standardized bias
 * the average confidence interval length
 * the average coverage of the CIs

Both these two methods are shortly discussed below.  

### IBMA
This approach is rooted into the meta-analysis literature. For every voxel, we calculate a weighted average using effect size estimates for the true population effect. There are several options for the weights, as well as calculating the CIs. Some of these are discussed in this paper.

### Third level GLM
As is standard in fMRI, statistical modeling in a single study consists of applying a two-stage procedure. First, individual time series are modeled within each subject using a General Linear Model (GLM) in each voxel. Then several subjects are combined in the second stage using again a GLM. The idea is to pull this through at a third stage in which studies are combined using a GLM. Intersubject, between-subject and between-study variability can be accounted for in each stage.

# Link with OSF
This project is hosted on the Open Science Framework and can be found using [this link][1].
 [1]: https://osf.io/t92bd/

# Structure of repository
I have not found time to come up with a decent structure.
The R files are listed in the main folder. A selection of reports can be found under _Reports_.
Some R Objects that contain analyses are listed under RObjects. These are merely used to speed up analyses.


# Getting started with Git
* https://git-scm.com/doc
* https://www.youtube.com/user/GitHubGuides/videos


# Contact
Webpage = [here][Han Bossier] or through [e-mail](mailto: Han.Bossier@Ugent.be).

[Han Bossier]: http://telefoonboek.ugent.be/nl/people/802001626303
[Freya Acar]: https://telefoonboek.ugent.be/nl/people/802001860820
[Ruth Seurinck]: http://telefoonboek.ugent.be/nl/people/801001629152
[Beatrijs Moerkerke]: http://telefoonboek.ugent.be/nl/people/801001453542


First Created: 12/01/16
