# Advanced topics in scRNA-seq analysis, May 2020


[ToC]
## Day 1

### Questions
#### Python and R
* **Q**: Can Python import sparse matrices via reticulate?
    * **A**: Yes - we will se this in the RStudio exercise (specifically, `Matrix::dgeMatrix` objects are converted by `reticulate` to SciPy CSC sparse Matrix objects)
* **Q**: How are classes handled
    * **A**: you need to write custom converter for non-standard class. (Not covered in exercise)
* **Q**: Is there a way to get help on a python function through R?
    *  **A**: In a python chunk: `help(function_name)`. From an R chunk, you can use the `py_help()` function.
* **Q**: What is the equivalent of %Rpush/%Rpull when using rpy2.robjects ? Are these magic functions always accessible?
    * **A**: %Rpush/%Rpull would still work here (after the `%load_ext rpy2.ipython`), and can be also applied to instances of `rpy2.robjects`.
* **Q**: Somewhat unrelated question, does running bash chunk inside RMarkdown requires any special set up or is it as smooth as python chunks as we saw?
    * **A**: Bash chunks (marked with {bash} are supported by RStudio/knitr, but not in an interactive session (e.g. running a line of bash using Ctrl+Enter). See e.g. here: https://bookdown.org/yihui/rmarkdown/language-engines.html

#### Quantification of gene expression
* **Q**: How does "selective alignment" compare to pseudo-aligners like kallisto?
    * **A**: It is similar in performance (high performance), but gives better accuracy (e.g. for complex cases that do require end-to-end alignment)
* **Q**: Is there a minimum UMI length for UMI deduplication?
    * **A**: The problem is UMI collisions, which depends on sequencing depth. 10x Chromium v3 chemistry has gone to 12bp UMIs, which makes collisions very unlikely.
* **Q**: Did you compare Alevin to other tools using simulated data?
    * **A**: Yes (not shown on slides, but in paper). In summary, Alevin works better also on simulated data, especially for genes with lower uniqueness.
* **Q**: How does alevin compare in terms of memory usage with other tools?
    * **A**: It uses less memory than e.g. STAR or CellRanger (even though it is faster). Memory usage has gone down in recent versions of Salmon (details in the paper).
* **Q**: Does the memory usage depend on the number of threads?
    * **A**: No, except for very small differences, it is largely independent.
* **Q**: Are there specific parameters to be used for single nucleus data?
    * **A**: Yes, it would probably make sense to include introns, which will be discussed by Charlotte (RNA-velocity). However, it does not seem to make a big difference in many cases for simple gene quantification (as opposed to RNA-velocity analysis)
* **Q**: How can Alevin be used for scATAC-seq data?
    * **A**: In principle, the genome can be indexed, but that will use a lot of memory (~20GB for a mammalian genome). In principle it should already work now, but this is something that is currently actively further developed.
* **Q**: Does Gencode have an advantage over Ensembl?
    * **A**: The underlying annotation is the same, but Gencode is preferred because it is self-consistent (transcript fasta files are consistent with the coordinates in the gtf file).
* **Q**: Does -2 always contain the R2 read?
    * **A**: No. Alvin expects -1 to always be the read containing the barcodes. If these are contained in the R2 read, you need to manually swap R1 and R2 (provide R2 to -1). You may also have to specify appropriate values to `--end`, `--umiLength` and `--barcodeLength`
* **Q**: Does `--noem` lead to multi-mapping reads to be dropped or distributed uniformly?
    * **A**: It would lead to dropping of multi-mapping reads, but only after UMI-deduplication, so it may still consider some level of multi-mapping and perform better than tools to completely ignore multi-mapping reads.
* **Q**: How precise should the expected cell number be? Is 1000 and 1200 any different?
    * **A**: Generally it will work even without giving that number. If the drop of barcode frequency ("knee-plot") is not very steep, it may however be difficult to auto-detect. `featureDump.txt` allows you to look at that distribution and choose an appropriate value. `--noQuant` can be used to just get that distribution without full quantification.
* **Q**: Does providing the genome as a decoy sequence during indexing improve performance?
    * Basically the idea is where best a read can align, apart from the transcriptome. If you provide whole genome as a decoy it does improves the performance of quantification. We have discussed in detail about this in the following paper. https://www.biorxiv.org/content/10.1101/657874v1. The SAF version of salmon here includes genome as a decoy sequence.

### Renku
* If you are from a swiss educational institution and you have the edu-id configured, you can log in with that.
* Interactive environment will be shutdown after 24h of inactivity
* Command to add course github: `git submodule add https://github.com/fmicompbio/adv_scrnaseq_2020.git adv_scrnaseq_2020`
* To avoid merge conflicts, consider the `adv_scrnaseq_2020` subfolder <span style="color:red"> read-only</span>
* To save any changes/material, remember to track them with `git add`
* To get latest updates: 
```shell=bash
cd adv_scrnaseq_2020
git pull origin master
```
Or:
```shell=bash
git submodule update --remote --merge adv_scrnaseq_2020
```
* unpushed changes found at `renku/autosave...` branch indicates that you had changes that were not addes/committed/pushed when you closed your environment, and were thus safed in the `renku/autosave/...` branch. If you want to recover these, you need to merge this branch into master (e.g. using git on the command line or the GitLab interface from renku in your browser)
* If you have a JupyterLab session running and instead want RStudio, replace `lab` with `rstudio` in the URL. Reverse is true when switching from RStudio to JupyterLab.

<img src="https://i.ibb.co/7yY36QT/image-2.png">

* In order to reduce memory usage, it may be useful to close e.g. RStudio (quit the session) when switching to JupyterLab

* If you want to run the same docker image locally, you can do that using:
    *  `docker run -ti --rm -p 8888:8888 -v ${PWD}:/work <image-name> jupyter lab --ip 0.0.0.0 --port 8888` (you can ignore the warning about the missing `/tmp/renku-env`).
    *  `<image-name>` is following the syntax: `registry.renkulab.io/<username>/<project-name>:<short-commit-sha>` where `<short-commit-sha>`  are the first 7 characters of the commit sha, e.g. `registry.renkulab.io/stadler.michael/adv_scrnaseq_2020:d5ee78a`.
    *  So putting things together, to run the same image that Rok had from a fork of the course project, you would do: `docker run --rm -ti -p 8888:8888 -v ${PWD}:/work registry.renkulab.io/rok.roskar/adv_scrnaseq_2020:d3ea747 jupyter lab --port 8888 --ip 0.0.0.0`, (assuming PWD is actually the checked out repository).
    *  This will require root (e.g. use `sudo docker ...`).
    *  You can see the list of images and tags here: `https://renkulab.io/gitlab/<user-name>/adv_scrnaseq_2020/container_registry`.
    *  The images will be available for download from renku for the forseeable future. Alternatively, you can clone the repo and run `docker build -t <image-name> .` to re-build the image locally.

### Combining R and Python

#### Reticulate in RStudio
* **Tip**: when using Rmarkdown, run the code top-to-bottom to render HTML-file to ensure output is correct.
* If you are prompted about whether you want to download and install Miniconda in renku, answer **no**.
* To enter the python REPL use command `repl_python()`. Then, to exit type `exit`.
* `reticulate` does not always pass by reference, copies can be created - take into consideration when working with large objects.
* you can use `py$x` from R to access Python variable `x`, or `r.y` from Python to access R variable `y`

#### rpy2 in JupyterLab
* **Tip**: `rpy2.robjects` (the high-level interface) is a good compromise between ease-of-use and efficiency

### Quantification of gene expression
* UMI length of ~12bp gives you a very low probability of collisions
* alignment vs. mapping: read mapping assigns reads to e.g. genes, but does not create a end-to-end base alignment (CIGAR string)
* UMI deduplication is non-trivial and one focus in the development of Alevin
* many tools drop ambiguous reads (multi-mapping reads), which are resolved by Alevin using an EM algorithm
* while genes with high uniqueness are also well handled by e.g. CellRanger, it systematically gets worse as uniqueness decreases, while Alevin shows good performance that is independent of gene uniqueness
* it may not be worth going to more than 12 or 16 threads, as Salmon may become IO-bound (limited by the speed of the storage system)
* **Inferential replicates**: Biological replicates are obtained by performing the experiment multiple times. Inferential replicates are multiple expression level estimates which Alevin can generate (`--numCellBootstraps` parameter, generates mean and variance estimates). It's an in silico approach to provide information on estimate uncertainty, by different assignment of ambiguous reads; how to use them in analysis is a topic of ongoing research (collaboration with Mike Love, imported into R using `tximeta`, visualization using `fishpond`, differential expression analysis using `Swish`). For example, you could ignore identified marker genes, if they have a very high inferential variance.
* Ongoing research: How to improve estimation by sharing information across cells:  https://www.biorxiv.org/content/10.1101/2020.04.10.035899v2



## Day 2

### Prepare your environment for the rna-velocity section:

* **Start an environment with 4Gb of memory**
* **Update the course material**
git submodule update --remote adv_scrnaseq_2020
* **Copy the Rmd file and figures to the working directory**
cp adv_scrnaseq_2020/rna-velocity/rna-velocity.Rmd .
cp -r adv_scrnaseq_2020/rna-velocity/rna-velocity-figures .


### Questions

* **Q**: Where does the $\gamma$ and $\beta$ parameters estimates come from
    * **A**: These are not known values, we estimate them from the data. Actually it is their ratio (the slope) that we estimate (assuming steady state is reached for some cells).
* **Q**: Going back to beta and gamma values, so with real data, we don’t actually need to estimate/set the exact values of beta and gamma, as long as we assume that some cells have reached equilibrium (top right corner), all we need is the ratio of gamma over beta to calculate velocities. Is that right?
    * **A**: Only under the assumption that some cells had reached steady state, yes (velocyto or scVelo-steady-state-model). Otherwise, they have to be estimated (scVelo-dynamical-model).
* **Q**: How do we handle genes that do not have introns or splice variants. 
    * **A**: Method is not applicable to genes without introns. Genes without splice variants are ok, as long as there is both exonic and intronic regions.
* **Q**: Experimental confirmation of profiles?
    * **A**: Validation is tough.  However results made sense in terms of the known biological mechanisms, e.g. velocity vectors were pointing from progenitor cells to differentiated cells when projected to a 2-dimensional embedding. RNA velocity was estimated to "look into the future" about 1-2 hours. 
* **Q**: Effect of single nuclei data on velocity estimation?
    * **A**: Open question but in summary the velocyto method was not explicitly developed for such data BUT seems to work also in that setting according to authors
* **Q**: I realize this analysis is typically done in highly curated genomes like mice or human, but I'm curious to know whether there were some estimations of the impact of the precision of the annotation on the overall results. For instance, by testing different versions of the genome or annotations.
    * **A**: Probably yes, method assumes knowledge of intron and exons, might perturb the system by "mismapping" and render erronous results. Recommend to inspect the results.
* **Q**: How many genes are concerned by those situations? roughly (overlapping genes)
    * **A**: If you throw away multimapping genes that could affect up 10-15% of the reads. But can also affect 40-50% of the genes, for well annotated genomes.
* **Q**: Are the full length technologies better for RNA velocity analysis? Having more info along the gene seems to be adventagus
    * **A**: You could imagine so but you also have _fewer_ cells, meaning it will be harder to estimate the actual trajectory.
* **Q**: in the same line as the question above, is it possible to do these kinds of analyses on non-umi non-droplet datasets like SMARTseq? 
    * **A**: It can be used, absolutely.
* **Q**: For the 3’ end data does the read length have any impact? 
    * **A**:  Yes it can, in general longer reads are better since they map better. Easier to see partial mapping (exon and introns) with longer reads.
* **Q**: from a practical point of view, I normally analyze my data on galaxy project, where alevin and starsolo are available. are there specific settings on galaxy that you need to take into account for velocity analysis?
    * **A**: Depending on which version that is installed `starsolo` should work, specify that you want the velocyto output. For `alevin` slightly more complicated, you need to fix the reference; if you have the option to provide your own fasta-file that would work but you would have to prepare an extended version. Creating the fasta-files and a Salmon index is covered in the RNA-velocity exercise. It is also described here is an Alevin tutorial: https://combine-lab.github.io/alevin-tutorial/2020/alevin-velocity/
* **Q**: Difference to "classical" pseudo-time ordering / trajectory estimation methods?
    * **A**: Those methods cannot get any information on the directionality of the trajectory. On the other hand the trajectories that can be estimated with velocity-based methods need to have dynamics that "live" on similar time-scales as the splicing dynamics.
* **Q**: What about cases where several orthogonal (or partially confounded) dynamic processes are taking place in the cell population?
    * **A**: If we have an idea of which gene-sets are involved in those different processes then we can use those to estimate distinct trajectories for each process. `scVelo` has a mode that tries to estimate the presence of distinct trajectories in **different** cell subpopulations. However, this assumes that those trajectories are not taking place within the **same** cell subpopulation.
* **Q**: Is there a QC measure that would tell you that velocity can be estimated?
    * **A**: Pre-filtering of genes with high-enough counts is built into velocity analysis tools. One could use similar approaches as for standard quantifications to get information on how reliable counts are.

 

### RNA velocity
#### Theoretical Session
* A way to estimate cell trajectories without having longitudinal data
* Notation:
    * $u(t)$ - unspliced
    * $s(t)$ - spliced
    * $\alpha$ - transcription
    * $\beta$ - splicing rate
    * $\gamma$ - degradation rate

Main system of equations:
\begin{equation} 
\begin{array}{ll}
& \frac{du(t)}{dt} = \alpha_k(t)-\beta u(t)\\
& \frac{ds(t)}{dt} = \beta u(t) - \gamma s(t) 
\end{array}
\end{equation}

* $\beta$ and $\gamma$ assumed to be constant over _time_
* $\alpha$ can change to the state that the cells are in
* The phase plot is a summary of the changes of the spliced and unspliced fractions of the RNA.
* Equlibrium means that the rates are zero (our system does not change) 
* Velocities ($\frac{ds(t)}{dt}$) are astimated as the vertical distance from the "steady state line" (slope $\gamma/\beta$) in the phase plot. 
* We model the system with relative time units (think of them as proportional to actual time)
* If we the system has not reached steady state, this could lead to overfitting unless accounted for in the modeling (advantage of `scVelo`)
* Correct definition of intronic regions is extremely important for robust velocity estimation.
* Different strategies of intron/exon read counting (e.g how are ambiguous regions handled) can have a strong impact on the final velocity estimates.
* Overlapping features or multimapping features can have large impact the inference, keep in mind that difference methods have different strategies for this. A lot of tools may just "_throw away_" reads overlapping such features. 
* Another difference is how different tools handle strandedness of genes (i.e whether reads coming from the opposite orientation as the annotated feature are counted). 
* With 3' data it's easier to say if they come from introns or exons than if they come from the spliced or unspliced version of a gene or transcript.


#### Exercise Session

* **Q**:What do you do with overlapping genes (on the same strand)?
    * **A**: We don't bother too much about them, but let the software take care of it (e.g., `alevin`) 

* **Q**: Why combine exons but not introns
    * **A**: you wouldn't gain anything from combining introns. But for exons it makes sense since genes may map to several of them.
    
* **Q**: How does `scVelo` compare to methods that do not rely on splicing information, for example `slingshot`.
    * **A** : These alternative methods tend to be based on distance metrics, rather than dynamical systems. In short, You embed your data in lower dimensional space and then try to construct a trajectory (one dimensional path, may be branched) along which your data points (cells) are arranged. These methods however do not really provide a notion of direction, start and end-points. But the trajectory based methods make fewer assumptions on the underlying system and are therefore to prefer in scenarios where the assumptions on time-scales do not hold, for example very rapid or slow changes in the system. One additional benefit with the methods built around dynamical systems and usage of exons/introns is that the phase portraits may also serve as "sanity checks", allowing us to assess how well the data fits the model.

* **Q**: Does Seurat contain a wrapper for RNA-velocity analysis?
    * **A**: Yes, it wraps the "old" velocyto (static model only, assuming steady-state), but not `scVelo` that also implements the dynamical model.
* **Q** : Where did the `.obsm["X_pca"]` and similar objects come from?
    * **A**: These were taken from the R object (have been pre-calculated). 
* **Q**: Can you comment on why 30 for `min_shared counts`? 
    * **A**: Always up for discussion, depends on the objective of your analysis. Hard to give a value that works for every case - also dependent on your data.
    
* **Q**: What if it's already normalized, and you enforce normalize again? 
    * **A**: it would not be the gene normalization method, but you would have to normalize the exon and intron matrices. If these are already normalized you can skip the normalization step.

* **Q** : When we say normalization here, are we referring to standard library-size normalization or is there more to it? 
    * **A**: We normalize the library size to the median.

* **Q**: How can we easily see genes that were filtered out? 
    * **A**: Just compare the genes that are left with all of the genes that were there before. Example by saving a variable with all the gene names in the original data before filtration and then comparing this to the current genes (after filtration)

* **Q**: It seems that the final interpretation heavily depends on how well the 2D visualization (tsne, umap, etc.) captures / distinguishes different populations? Is that right? For example if early and mid round populations are not well separated, then velocity vectors will be circular or something? 
    * **A**:  Yes, this is true. In theory, the most interpretable space to visualize these are in PCA space, where we have a proper interpretation of distances. Be very careful with tSNE and UMAP; good idea to assess the results with several methods.

* **Q**: How did you interpret the red cluster in the middle? 
    * **A**: (from Sina) Regarding the red isolated cluster, this paper suggests it could be an artifact of distortion in dimensionality reduction: https://www.biorxiv.org/content/10.1101/689851v3 

* **Q**: Can you import these vector coordinates from python back into r and visualize them easily there? 
    * **A**: yes that is possible, the vectors are accessible from the `AnnData` object. 
    
* **Q**: When a gene is said to be "downregulated" what does this mean?
    * **A**: it means that the abundance will decrease over time, it's a consequence of degradation being higher than production. This takes us from the steady state.

* **Q**: Can we examine genes that have not been fitted?
    * **A**: No not really, but we can visualize their expression.

* **Q**: Would it be possible to do RNA-velocity analysis in a system with oscillating genes?
    * **A**: In principle yes. The cells would be expected to go around the phase plot multiple times. Pseudotime estimation might be difficult (`scVelo` might not know which cycle a cell is in).

* **Q**: Would a trajectory inference result roughly agree with scVelo?
    * **A**: Yes, in this particular case it likely would, although the inferred trajectory (depending on the tool used) would not have directionality information.

* **Q**: How useful is it to try to do velocity analysis on datasets where you don’t expect clear trajectories? 
    * **A**: In that case it may not really give you much. For example, when applied to dataset with just differentiated cell types, you get disconnected clusters which do not really say much.


* **Q**: If we have 3 single cell experiments, do you suggest to run it separately or should we run it on the complete data set simutaneously?
    * **A**: It will not really take a normalized integrated data set, hence if you need to regress out strong batch effects you might therefore run into problems. No "established" truth or best practice. If you do not expect batch effects then try using all in one run. On the brink of what we currently know.


<hr>

* **Disclaimer** : this is not "the only" approach - area of active research.
* `scVelo` now offers the ability to fit multiple models (and evaluate whether this approach is to prefer) to certain cell subsets.
* Usage of the EM algorithm in `alevin` may result in non-integer values. Important to take into consideration when designing your worklflow.
* Number of genes that you select (`n_top_genes` in `filter_genes_dispersion` function) should to some extent relate to the complexity of your tissue. 
* The **Dynamical** approach has a tendency to increase run time quite significantly, given how it solves the system of equations. Often gives you _better_ and more _coherent_ results, usualy worthwhile the extra time.
* Solving the equations is not that easy since we do not know which timepoint each cell corresponds to. If we knew this information we could simply plug in these values and then solve the equations. To circumvent this issue an EM algorithm is used - iterative inference.
* $\alpha,\beta,\gamma$ are usually intialized with values from the steady state model.
* Trajectory inference of a gene is independent of the other genes. 
* To obtain the angle between two vectors $x$ and $y$ we can simply use the fact that:
\begin{equation}
    x\cdot y = ||x||||y||\cos(\theta)
\end{equation}
Where $\theta$ is the angle between the two vectors. This can then be used to get the angle between the displacement and velocity vector.

* For low dimensional representation of the velocities we convert the cosine similarities to transition probabilities ($\tilde{\pi}_{ij}$). These can be interpreted as the probability of a cell ($i$) transitioning into the state of another cell ($j$). Consider the low dimensional representations as a "weighted average" of the displacement vectors in the low-dimensional embedding, weighted by the probabilities of moving towards the other cells.
* the streamline plots aims to represent the underlying field of the phase space.

### Spatial transcriptomics

* focus on Visium technology (10x-genomics commercialized version of the 2016-"Spatial Transcriptomics" method)
* hexagonal grid of ~5000 spots in ~6.5mm-by-6.5mm square (resolution ~100um, 1-10 cells), polyT capture oligos with spatial barcode (millions of identical oligos per spot)
* data is processed using SpaceRanger (similar to CellRanger)
* Visualization: Use DimRed techniques to move to 3D then convert to RGB space by mapping each dimension to a color channel.
* Preprocessing, filtering and normalization similar to the approaches used for SCdata (filter spurious spots, non-detected genes, normalize for spot library size, remove mt/rib genes). Might need to include chip as covariate as they tend to exhibit fluctuating patterns.
* spots measure 1-10 cells and thus can be mixtures of cell types -> spot clusters are not "pure" cell types; spot clusters can be thought of as variable-frequency mixtures of pure cell types
* Spatial mapping of cell types by integration of transcriptomics data (`stereoscope` tool): https://www.biorxiv.org/content/10.1101/2019.12.13.874495v1. The implementation accounts for asymmetries (cell types only present in one data type, either spatial or single-cell data). One advantage is that the approach is independent of specific marker genes.
* Often useful to filter out mitochondrial and ribosomal genes if we're not specifically interested in them.
* Beware that the way you transform and scale your data strongly affects how it will look when visualized. 
* publicly available visium data: So for large sets of Visium data there are not all that many resources. However, data will soon be published here: https://hdca-sweden.scilifelab.se/tissues-overview/ and for the old ST technique (1k spots) you can find sets with more replicates from more individuals here: https://www.spatialresearch.org/resources-published-datasets/




### Questions

* **Q**: Is it possible to use cell type deconvolution techniques to infer the composition of a spot?
    * **A**: Yes - we will see how to use single-cell data to make inferences about the spatial data. 
* **Q**: Does it make sense to use the spatial information as part of the input for the clustering? 
    * **A**: Can be done, depends on what you want to get out. Not used very often - risk that the spatial information drives the clustering too much.
* **Q**: How would you "do an affine transformation to unit cube" for visualization?
    * **A**: After embedding in 3D, divide each coordinate by the max value, to get all three coordinates into [0,1] (so they can be used as RGB color channels)
* **Q**: How do you address smoothing gene expression over 2D?
    * **A**: Gaussian process or loess (for expression as a function of distance)
* **Q**: Are there methods for spatial gene set testing?
    * **A**: Estimate an enrichment score for each gene set in each spot.
* **Q**: Is spatial data more or less sparse than 10x gene expression data?
    * **A**: Typically a bit more sparse - more things that can disturb the capture of the genes. 
* **Q**: Origin of variations between different chips?
    * **A**: Origin of variation is a bit unclear. But you see the separation by chip ID clearly in e.g. UMAP. Each chip corresponds to a different section - adjacent sections can often be aligned, but not if sections are taken very far apart in the tissue sample.
* **Q**: Do you know if any of these spatial techniques have been combined with expansion microscopy already?
    * **A**: Yes, there was a recent preprint: https://www.biorxiv.org/content/10.1101/2020.05.13.094268v1
* **Q**: When looking at expression as a function of the distance from a cluster, what about the cells within a cluster?
    * **A**: They are excluded in this analysis. One could also imagine considering only these cells, or multiplying their distance to the boundary of the cluster by -1.
* **Q**: What R packages are there for spatial transcriptomics?
    * **A**: Seurat, STUtility (https://ludvigla.github.io/STUtility_web_site/, wrapper around Seurat)
* **Q**: Can one apply velocity to spatial data?
    * **A**: This is discussed e.g. here (for MERFISH data): https://jef.works/blog/2020/01/14/rna_velocity_analysis_tutorial_tips/ (paper: https://www.pnas.org/content/116/39/19490). The idea is that the nuclear mRNAs = unspliced, and cytoplasmic mRNAs = spliced. One issue for ST/Visium data is the mixture of different cells in a spot, and that they may not all be in the same cell state, which makes the inference difficult since only the joint expression profiles are observed.
* **Q**: Does the interpolation (in notebook 3) use also the HE data?
    * **A**: In this case no, but it could be used. 
* **Q**: Experience with nanostring spatial data?
    * **A**: Seems to be gaining traction.

## Day 3
### On-disk data

**Update the teaching material**

```shell=bash
git submodule update --remote --merge adv_scrnaseq_2020
```
**copy the on-disk-data folder to our working directory**

```shell=bash
cp -R adv_scrnaseq_2020/on-disk-data .
```

`renku` specs:
* Environment : `rstudio`
* #CPU : $\geq 1$
* Memory : 8GB
* LFS-data box : Checked (**Important**)

<hr>

#### Questions

* **Q**: How do we handle empty columns
    * **A**: It will just repeat the column pointer. Here is an example in R:
    ```
    library(Matrix)
    set.seed(2)
    x <- sparseMatrix(i=sample(4, 2), j=sample(4,2), x=sample(100,2), dims=c(4,4))
    x # columns 1 and 3 are empty
    x@i # row indices
    x@p # column pointers (note the repetitions)
    x@x # non-zero elements
    ```
* **Q** : Could you also use linear indices to remove one of the two indices? Or even differences between these linear indices.
    * **A**: You could in theory do that, but it not the convention within the community.

* **Q**: In the examples with chunks as a whole column, is that bad when accessing the diagonal of the matrix.
    * **A**: Yes, it would be a relatively bad strategy since you would have to read the whole matrix.
    
* **Q**: How does the chunk strategy apply to _sparse matrices_
    * **A**: TBC
* **Q** : are there any standards regarding what files you use as input, put into you HDF5 files.
    * **A**: you do not usually put other files in the HDF5 files. Some formats have more of a defined/enforced structure and it's specific character.
* **Q**: Any idea why the loom-file requires a dense matrix.
    * reading a dense matrix is an easier task than reading a sparse matrix; might be due to that.

<hr>

#### Lecture Notes

* Memory requirements scale linearly for dense matrix representations
* Sparse matrix representation utilize three vectors to represent the data
    * Values : non-zero values
    * Row Indices : which element is _value_ located
    * Column Pointers : at which element in _row indices_ new column starts
* In HDF5 every element - even top - can all hold metadata $\rightarrow$ "self-describing" data.
* HDF5 do not store data sets as contigous data but rather in chunks. The purpose being to avoid the issue of efficient reading along one dimension but _very_ poor efficiency along the other.
* Chunks in HDF5 not not necessarily have to be square. Picking the dimensions is of importance when designing for high efficiency.
* There is a tradeoff : entire cunk needs to be read (restrict size), but there is also an overhead to finding chunks (restrict number of them). Can have an huge impact on the performance if chosen correctly (or poorly)
* Recommend to not really use other filters, given how there is no guarantee that users have the necessary backends.
* HDF5 "chunk-strategy" allows us to use dense matrices, rather than sparse.
* Effectiveness of compression $\propto$ chunk size
* The ideas may be transferred to sparse matrix representations as well. Each vector can be chunked and compressed independently. Not quite as much gain as for dense matrices but still offers some benefits.
* Accessing data from disk is **much** slower than access from memory $\rightarrow$ negative impact on processing time.

<hr>

### Exercise session

#### Questions:

* **Q**: Can we `head` and hdf5 file. 
    * **A**: No not really i the same way
* **Q** : can we store images?
    * **A**: yes, it should be possible.
* **Q** : Error
```r
Error: Error in h5checktype(). Argument not of class H5IdComponent. 
h5writeAttribute(gene_names, "../data/on-disk-data/my_hdf5.h5", name="gene_names") 
```
* **A**: Thanks for noticing this, it's gone on my [list ](https://github.com/grimbough/rhdf5/issues/62) of things to update in rhdf5.  Your usage makes sense and it's the software that's at fault.  The error is because `h5writeAttributes()` works slightly differently and expects to be given an object that represents an open HDF5 dataset rather than just a path - we didn't cover this in the lecture; HDF5 has a lot levels of functionality.

You can achieve this at the moment this by attributes to an object in R before writing it to the HDF5 file e.g.

```r
> ex <- matrix(1:10, ncol = 2)
> attributes(ex) <- list(scale = "minutes")
> attributes(ex)
$scale
[1] "minutes"
> h5write(ex, file = "example_file.h5", name = "mat_with_att", write.attributes = TRUE)
> h5readAttributes('example_file.h5', name = "mat_with_att")
$scale
[1] "minutes"
```


* **Q**: If you want to add a new group would you have to restart then?
    * **A**: You can have groups within groups, but for the dataset can not have groups.
* **Q**: We only operate with chunks when we use matrices - is that true?
    * **A**: No, not necessarily. HDF5 can store files with up to 16 dimensions. Datatype may always vary, e.g., strings.
* **Q**: You mention that h5 support `gzip` natively, could you exemplify some more.
    * **A**: it has `gzip` and one that the h5 group maintains `szip`, plus a few others (there's a concise summary [here](https://portal.hdfgroup.org/pages/viewpage.action?pageId=48808274)) . They also support 'plugins' that are not distributed with the main HDF5 library but can be used dynamically when needed - useful when you have very intrinsic knowledge of you data.  There's some more details from HDF5 group [here](https://portal.hdfgroup.org/display/support/HDF5+Filter+Plugins) and from Mike in the [rhdf5filters](https://github.com/grimbough/rhdf5filters) package.
* **Q**: Is there any way to directly write metadata (e.g. a new meta data column) to an hdf5 file? 
    * **A**: You can either write a new dataset containg this meta-data and give it a sensible name that describes what it is. Or you could write an attribute to the dataset that the metadata relates to. 
* **Q** : the time to read depends not only on the amount of columns but also on the position where these columns are?
    * **A**: Not expected behaviour, should not be dependent on position - as long as the indexing is continous.Picking random cells is suboptimal

* **Q**:Maybe I missed it but how the chunk size is encoded/choose when you create a file? 
    * **A**: We did not specify that, will look at it later (towards end of 03-on-disk.Rmd). The chunk size in our set is 100 (prior knowledge)    
    * The HDF5Array package will let you verify the chunk size of an existing dataset.  We didn't really introduce HDF5Array, but the code below will return the chunk dimensions: 
    ```
    > HDF5Array::chunkdim( 
            HDF5Array(file = "data/on-disk-data/brain100k.h5", 
                      name = "/counts_matrix" ) 
        )
    [1] 100 100
    ```
    * It is also possible to use the rhdf5 package, but requires a lot of functionality we didn't explore.    
    * Setting the chunk size can be done when creating a dataset with `h5createDataset()` and specifying the `chunk` argument.  Chunk size is immutable once a dataset has been created.

* **Q**: How does HDF5 handle `NA` values?

    * **A**: This is implementation specific as the HDF5 specification doesn't have a concept of NA.  This is an outstanding issue in rhdf5 (e.g. [#58](https://github.com/grimbough/rhdf5/issues/58), [#61](https://github.com/grimbough/rhdf5/issues/61)).  `NA` strings are handled OK, but not so well for integer or numeric data types.

### DGNs

#### Preparation
**Update the teaching material**
```shell=bash
##!!Only needed if you have not updated today**
git submodule update --remote --merge adv_scrnaseq_2020
```
**copy the on-disk-data folder to our working directory**

```shell=bash
cp adv_scrnaseq_2020/DGNs/Keras_example.Rmd .
cp adv_scrnaseq_2020/DGNs/DGNs_exercise.Rmd .
```

#### Lecture notes:
* Tensorflow is the most widely used computational back-end for machine learning
* Tensors are (n-dimensional) arrays, layers are computations (transformation functions) applied to tensors;
  combined, they form the computation graph
* Keras is a higher-level API, providing convenience wrappers for commonly used layers or graphs (R and Python version available)
* Deep Learning:
    * Model that takes input and transforms it to an output via successive layers of increasingly abstrand/meaningful representations
    * "meaningful" depends on the task at hand and is "learned" by the model during training
    * "deep" refers to the (many) layered representation
* Training:
    * important: definition of loss function (how to measure model performance)
    * iteratively update model parameters such the the loss function is optimized (e.g. the error minimized) -> done by optimizer
    * the error observed at the model output is back-propagated through the model towards the input, and parameters are updated based on the local gradient
    * the training data is processed in chunks (*batches*) and weights are updated after each batch; one loop through the training data is called an *epoch*; many epochs of training are usually performed until convergence
    * if batch size is too small, model performance may fluctuate a lot; if it is too large, it may not fit into e.g. GPU memory; learning rate and batch-size parameters are associated and have to be defined appropiately
* Deep learning revolution was spurred by 1) hardware improvements (e.g. GPUs), 2) software improvements (algorithms) and 3) availability of large high-quality datasets
* Autoencoders
    * unsupervised (need no labelled data)
    * "hourglass" topology (encoder - latent code - decoder) with "bottleneck" in the middle (latend code)
    * bottleneck forces model to extract only the most salient features of the input, that are most needed to reconstruct it (compression)
    * variational autoencoder: learn latent code distributions (characterized by means and variances) instead of point estimates -> allow to interpolate and explore unobserved samples
* Generative Adversarial Networks (GANs)
    * Generator (creates fake samples) and Discriminator (tries to classify if input is fake or real) components
    * objective is to create fake data that cannot be discriminated from real data
    * hard to train, notoriously unstable (mode collapse), but can create very realistic fake datasets, and are great to do inference, e.g. prediction of perturbation effects (generative network)

#### Questions:
* **Q**: How does relu handle negative values?
    * **A**: relu is defined as: $y = max(0, x)$; when a negative output is provided as input to relu it will produce a zero value - which is not equivalent to "negative values" not being compatible with the activation function. This behavior to some extent mimics the behavior of real neurons where a certain voltage is required to propagate the action potential (here the threshold would be value > 0). relu is not the only activation function, other alternatives exists e.,g Leaky relu, elu and tanh.
* **Q**: What is a Dense layer (Keras example on slides)?
    * **A**: There are many possible layers pre-defined in Keras. The dense layer is a layer that connects all nodes from the current layer with all nodes from previous and next layers. See these links for available Keras layers: [in Python](https://keras.io/api/layers/), [in R](https://keras.rstudio.com/articles/about_keras_layers.html)
* **Q**: How are model weight updated during training for 1 epoch?
    * **A**: Weights are update after each batch, so depending on batch-size, multiple times per epoch.
* **Q**: How do you choose hyper-parameters of your model (learning rate, number of layers/nodes per layer, dropout rate, ...)?
    * **A**: No easy recipe - learn from experience (your own and others, e.g. from literature/internet). A lot is also trial and error. Tools that sweep hyper-parameter-space to find optimal values are hard to use (too computationally expensive)
* **Q**: How does the dropout layer make training more robust / overtraining less likely?
    * **A**: A dropout layer sets a random subset of weights to zero - by that, it forces the network to have built-in redundancy and does not allow it to focus on non-general features in your input (e.g. dust or scratches on some images). Other techniques are available to avoid overtraining, e.g. L1 or L2 regularization.
* **Q**: Can you still use the same training (optimizer) for variational autoencoders as for standard autoencoders?
    * **A**: Yes.
* **Q**: Can Keras read data from disk?
    * **A**: Yes, including HDF5 - however I have not used it myself.
* **Q**: Would you suggest to do imputation usind autoencoders before scRNA-seq data processing, e.g. for cell type labelling?
    * **A**: This should work and make labelling easier, assuming that your model is trained well.
* **Q**: The autoencoder in the exercise achieves a reconstruction accuracy of ~0.5 (correlation) - would we not expect better?
    * **A**: Not necessarily - remember that the reconstructed data does not have the dropouts and the high-magnitude outliers anymore, and thus does not really look exactly like the input anymore.
* **Q**: When doing batch correction using latent-space-arithmentics, could you do it if your batches do not contain the same cells?
    * **A**: There are limits - the latent-code delta will also capture non-technical sources of variance that discriminate your badges, e.g. differences in cell type composition. If you have a common cell type, you could estimate the delta vector only on that and apply it to all cells. There are other more sophisticated ways to do batch correction using variational autoencoders (beyond the scope of this exercise).
* **Q**: Could you use an autoencoder for dimensionality reduction?
    * **A**: Yes, this has actually already been done (see also references given in slides and exercises).
* **Q**: Can you train an autoencoder model on a very small/limited dataset?
    * **A**: It may work by first training the model on larger datasets that may be similar to your target data, to get a pre-trained model that you then refine on your small target dataset.
* **Q**: Could you pretrain on 10x and refine on C1 data?
    * **A**: No, probably that will not work well because of the large amount of technical differences between these two technologies.