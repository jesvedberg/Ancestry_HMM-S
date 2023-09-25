# Ancestry_HMM-S Manual  

## Jesper Svedberg

### Quick overview

Ancestry_HMM-S (AHMMS) is a program designed to infer adaptive introgression from population genomic data. This README.md file contains a short user manual. The paper describing this software and our work to validate it has now been published in *Molecular Biology and Evolution* and can be found here: https://academic.oup.com/mbe/advance-article/doi/10.1093/molbev/msab014/6120794 . AHMMS is based on the software Ancestry_HMM, and more information about input file types etc. can be found at: https://github.com/russcd/Ancestry_HMM

### Installation

AHMMS can be installed either through Bioconda, or by compiling it from source code.

#### Bioconda

[Bioconda](https://bioconda.github.io/) is a repository of bioinformatics software for the [Conda](https://docs.conda.io/en/latest/) package manager. If you have [Conda installed and Bioconda set up properly](https://bioconda.github.io/user/install.html), you can install AHMMS using the following command:

        $ conda install ancestry_hmm-s

Bioconda is the easiest way to install AHMMS, especially if you do not have root access.

#### Download and compile

The easiest way to build AHMMS from source is to clone this Github repository and compile it:

        $ git clone https://github.com/jesvedberg/Ancestry_HMM-S.git
        $ cd Ancestry_HMM-S/src/
        $ make

#### Dependencies

If building AHMMS from source code, make sure the C++ linear algebra library armadillo is installed. More information and detailed download instructions can be found here, http://arma.sourceforge.net/. Armadillo can also be installed on OSX using homebrew

        $ brew install armadillo

On Ubuntu using apt-get:

        $ sudo apt-get install libarmadillo-dev

If you do not have root access, we recommend installing it using conda/miniconda: https://docs.conda.io/en/latest/

        $ conda install -c conda-forge armadillo

It is also recommended that users install the google-perftools package as compilation using tcmalloc tends to decrease runtimes, sometimes substantially. However, this is not necessary to use the software.

google-perftools can be installed using homebrew, apt or conda, in a similar way as described above. Conda example:

        $ conda install -c conda-forge gperftools

More information on compiling AHMMS using libraries installed using conda can be found in `src/Makefile`.


### Usage

For a list of options and arguments, see the built in help message:

        $ ahmm-s --help

This will print the following:

        ahmm-s usage:

                required:
                        -i [string]
                                input file name
                        -s [string]
                                sample id and ploidy file
                        -p [int] [int] [float]
                                ancestry pulse with format, ancestral population, time,
                                and proportion of final ancestry from this pulse
                        --ne [int]
                                effective population size of the admixed population

                select one of the following working modes:
                        --gss [int] [int] [int] [float] [float]
                                golden section search for optimal selection coeffient at each site.
                                parameters: chromosomal position start, stop, step, selection coefficient start, stop
                        --grid [int] [int] [int] [float] [float] [float]
                                calculate likelihood ratios in a grid.
                                parameters: chromosomal position start, stop, step, selection coefficient start, stop, step.
                        --site [int] [float]
                                calculate likelihood ratios for a single value of s at a single site.
                                parameters: chromosomal position, selective coeffient

                optional:
                        --help
                                print this help statement
                        -g
                                samples are specified with genotypes rather than read counts
                        --chr [string]
                                specify chromosome that will be analyzed
                                (only necessary when there are multiple chromosome in input file)
                        --chr_win [int] [int]
                                limit region on chromosome that will be analyzed
                        --gss_precision [float]
                                specify precision in finding optimal value of s using golden section search. default: 1e-5
                        --unit_coords
                                unit for start and stop position in grid and gss search can be defined as chromosome
                                coordinates rather than line in file. default off
                        --window [string] [float]
                                specify size of Markov chain in percent or Morgans.
                                "p 10" extends the markov chain 10% of chromosome length on each side of selected site.
                                "m 0.1" extends the windows 0.1 Morgan on each side of the selected site.
                        --traj [int]
                                change algorithm for generating selection trajectories.
                                4: 4-point approximation, 3: 3-point approximation (legacy option, not recommended),
                                default: forward iteration.
                        --stochastic
                                enables the stochastic method for generation selection trajectory.
                                (Experimental. Slow. Use for small values of s.)
                        --stochastic_reps [int]
                                specifies number of simulations for the stochastic trajectory algorithm.
                                default: 10000
                        --full_selection_space
                                turns off optimization of the selection coeffient search space. (Experimental)

### Usage concepts

AHMMS is a program designed to infer adaptive introgression from population genomic data. As input the following files and parameters need to be specified:

* A data file containing genomic data from a population (`-i filename`)
* A ploidy file specifying the ploidy of all individuals in the datafile (`-s filename`)
* Population size (`--ne int`)
* The time of the introgression event in generations, and the introgression fraction as a decimal number (see below for format).
* The analysis mode, with the choice of Golden Section Search, Grid Search and Single Site mode.

Further optional parameters for controlling the software are available as seen in the help message above, but default values are in most cases recommended. We do recommend using the following non-default parameters for improved performance, most likely without a loss of precision:

        --traj 4
        --window p 10

The `--traj 4` flag changes the method for calculating expected transition rates to a fast and accurate 4-point approximative method and the `--window p 10` flag changes the length of the Markov chain to only extend 10% of the chromosome length in each direction going away from the focal site, instead of including the full chromosome. Doing this will speed up computation time and is generally not expected to lower the accuracy of the method, but depending on for instance the density of variable sites in your data, you may want to make the window larger. Another option is to specify the window size in Morgans instead of as a percentage of the chromosome length. You can for instance set a window size of 2*0.1 Morgan using `--window m 0.1` (0.1 Morgan on each side of the focal site).

### Specifying time and size of introgression

Specifying the time since introgression and the size of the introgression pulse is done with the `-p` flag. AHMMS assumes that there is a single introgressive pulse that happened once in the past, from a donor population into a receiving population. These two populations are specified by using the `-p` flag twice. The `-p` flag needs the following parameters:

        -p POPULATION_ID TIME INTROGRESSION_FRACTION

`POPULATION_ID` is set to either `0` for the donor population, or `1` for the receiving population. `TIME` specifies the time since introgression in generations for the donor population. For the receiving population you can specify an arbitrarily large number (100000 is a good choice) to indicate that it is the receiving population that has been present for an arbitrarily long time. (The reason for this somewhat convoluted way of specifying time is that AHMMS is based on Ancestry_HMM which can handle more complex introgressive scenarios). `INTROGRESSION_FRACTION` specifies the size of the introgressive pulse. For the donor population it should be specified as the size as a fraction of 1, and for the receiving population it should be specified as one minus the fraction. If you for instance want specify a 10% introgressive pulse that happened 100 generations ago, you use the following command:

        -p 1 100000 0.9 -p 0 100 0.1

**WARNING! Bug in interpretation of `-p` parameter**

We have received user feedback showing that there is a bug in how the `-p` parameter is interpreted. Unless you specify population 1 before population 0, you will get nonsense values. In other words, please specify your population in the following way:

        -p 1 100000 0.9 -p 0 100 0.1

and not

        -p 0 100 0.1 -p 1 100000 0.9

We will fix this bug in a future update of Ancestry_HMM-S, but for now the program will work properly as long as you specify `-p 1` before `-p 0`.

### Analysis mode

AHMMS has three different analysis modes: Golden Section Search, Grid Search, and Single Site Mode.

**Golden Section Search** will find the selective coeffient at each site included in the analysis which produces the highest likelihood ratio. It does this through a simple hill climbing algorithm (golden section search). This is the recommended mode for inferring adaptive introgression from real data. To activate Golden section search, use the `--gss` flag with the following parameters:

        --gss WINDOW_START WINDOW_END STEP_BETWEEN_SITES MINIMUM_S MAXIUMUM_S

With WINDOW_START and WINDOW_END you specify the start and end point of the sites that you want to analyze. This is different from the `--chr_win` flag, which specifies the start and end position of the sites that will be included in the Markov chain. In order to include all sites in a dataset, set the start to 1 and the end to number of sites-1. You can also specify start and end in term of chromosomal coordinates (in bp). To do this you have to also use the `--unit_coords` flag. STEP_BETWEEN_SITES specifies the number of steps between all sites in the dataset that will be analyzed. This can be useful if you want to speed up your computation time, but it comes at a cost of precision in identifying the position of candidate sites. To analyze all sites set it to 1 and to for instance analyze every 10 sites, set it to 10.

With MINIMUM_S and MAXIUMUM_S, you specify the search space for inferring the selective coeffient s. We recommend setting these values th 0.001 and 0.15 as a start. In actuality, AHMMS will automatically cap the maximum value to whatever value of s that is expected to cause the site to reach a frequency of 0.99 in the specified time since introgression. This cap is used to prevent division-by-0 bugs if a site goes to fixation. You can turn off the cap with `--full_selection_space`, but be warned that this can cause AHMMS to report the likelihood ratio as `nan` under certain outlier conditions. Also, even using this flag, AHMMS will generally not report selective coefficients that are above the standard cap value.

Golden section search will print the selective coeffient that generates the highest likelihood ratio, together with the likelihood ratio for each analysed site to STDOUT. Column order: site, inferred selective coefficient, likelihod ratio:

        64154   0.00102039      -0.0422059912235
        71641   0.00102038967342        -0.0365331203211
        77484   0.00102038967342        -0.0374926163349
        ...

**Grid search** will analyze a set combination of sites and selective coeffients in a grid. In requires the following parameters:

        --grid WINDOW_START WINDOW_END STEP_BETWEEN_SITES MINIMUM_S MAXIUMUM_S STEP_S

This works in the same way as for Golden section search, with the difference you now specify the step length for all values of s to be tested. For instance

        --grid 1000 2000 10 0.001 0.01 0.001

will at every 10 sites between the 1000:nd and 2000:nd variable site calculate the likelihood ratio for values of s between 0.001 and 0.1 with step length 0.001 (0.001, 0.002, 0.003 ... 0.098, 0.099). The output will be printed to STDOUT with the same columns as described above:

        143625  0.001   0.0623890906572
        143625  0.002   0.113553284202
        143625  0.003   0.156187036075
        143625  0.004   0.18824519543
        143625  0.005   0.210872191936
        143625  0.006   0.21934426995
        143625  0.007   0.214521315414
        143625  0.008   0.194907535799
        143625  0.009   0.157839748077
        149168  0.001   0.0764014790766
        149168  0.002   0.141947525088
        149168  0.003   0.199271192774
        149168  0.004   0.246083804406
        149168  0.005   0.284079186153
        149168  0.006   0.307598042302
        149168  0.007   0.318045564461
        149168  0.008   0.314002908301
        149168  0.009   0.292439335957
        ...

Grid search is useful for visualizing the likelihood surface of your data (see figure 1D in our paper for an example), and can help when exploring your data, but it is generally slower than Golden section search for actually identifying the optimal values of s.

**Single site mode** will output the likelihood ratio for a single combination of a site and a selective coeffient. Usage:

        --site SITE S_COEFF


### Usage example

We have included a simulated dataset that can be used for testing if your installation of AHMMS is working properly. in the `example/` subdirectory, you will find the files `example.data` (containing genotype data), `example.ploidy` (containing ploidy information) and `readme` (usage information). The data is from a simulation of an introgressive scenario where a 1% introgression pulse took place 200 generations ago. At position 5,000,000, a locus with a selection coeffient of s=0.05 is located. To run this data use the following command:

        $ ahmm-s -i example.data -s example.ploidy  -p 1 1000000 .99 -p 0 200 0.01 --ne 100000 --window p 10 --traj 4 --gss 1 20360 10 0.001 0.15 > example.gss_output.txt

This will run Golden section search at every 10 sites across the chromosome, using the 4-point approximation method, with a window size of 10% of the chromosome and save the output to the file `example.gss_output.txt`. To visualize the output, please use your favourite plotting software.


### Preparing data to be used with AHMMS

AHMMS uses the same file format for genotype data as Ancestry_HMM. See https://github.com/russcd/Ancestry_HMM for further details.

To infer adaptive introgression you need the following data:

* Genotype data from
  * An admixed focal population
  * Two unadmixed parental populations
* A recombination map (optional)
* Simulations of a neutral introgression scenario to determine determine a likelihood ratio cutoff.

AHMMS has been validated using low coverage pileup data, but it is also possible to use genotype data, or data from a pooled sequencing experiment. Including a detailed recombination map will improve the power to detect adaptive introgression, but if one is lacking, a reasonable flat recombination rate may be enough to detect strong outliers. A simple Python script for converting VCF files to this format is found at `scripts/vcf2ahmm.py`. A short user manual for the script can also be found at `scripts/readme.md`.

Some parameters that are required by AHMMS has to first be estimated using other software. We recommend estimating the time since introgression and the size of the introgression pulse using Ancestry_HMM, though other similar software may also work well. You will also need to specify a population size, but our validation work has shown that AHMMS is not particularly sensitive to misspecification of this parameter, and a reasonable guess may be good enough.

See Schumer et al. (2020) https://onlinelibrary.wiley.com/doi/abs/10.1111/1755-0998.13175 for a detailed description of data requirements and a pipeline designed for performing simulations and inferring introgression using Ancestry_HMM. Much of this information is transferrable for using AHMMS as well. Also see Corbett-Detig & Nielsen (2017) and Medina et al. (2018) for more details on Ancestry_HMM.

### References

Svedberg, J., Shchur, V., Reinman, S., Nielsen, R., and Corbett-Detig, R. (2021). Inferring Adaptive Introgression Using Hidden Markov Models. Molecular Biology and Evolution. https://doi.org/10.1093/molbev/msab014

Further references:

* Corbett-Detig, R., and Nielsen, R. (2017). A Hidden Markov Model Approach for Simultaneously Estimating Local Ancestry and Admixture Time Using Next Generation Sequence Data in Samples of Arbitrary Ploidy. PLoS Genet 13.
* Medina, P., Thornlow, B., Nielsen, R., and Corbett-Detig, R. (2018). Estimating the Timing of Multiple Admixture Pulses During Local Ancestry Inference. Genetics 210, 1089–1107.
* Schumer, M., Powell, D.L., and Corbett‐Detig, R. (2020). Versatile simulations of admixture and accurate local ancestry inference with mixnmatch and ancestryinfer. Molecular Ecology Resources n/a.
