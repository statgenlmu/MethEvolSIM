
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MethEvolRSIM

<!-- badges: start -->
<!-- badges: end -->

Software for simulation of DNA methylation dynamics on different genomic
structures under the basic model for the evolution of DNA methylation
along genealogies. The methylation state of each CpG dinucleotide is
represented as 1 for unmethylated, 2 for partially methylated and 3 for
methylated.

# Installation

You can use the development version of MethEvolRSIM from git. The
development is being implemented with the `devtools` package. It can be
loaded and checked using the following:

``` r
library(devtools)
#> Loading required package: usethis
load_all()
#> i Loading MethEvolRSIM
```

`check()` runs the unit tests. The unit tests that test the functions
used for writing output files: - *validationResults.txt*: Is generated
in case an error relative to the characteristics of the objects happens
and contains both the error and the object.

# Manual of Use

## Structural distribution of CpG sites

The initial spatial distribution of CpG sites is given by the user. A
data frame with 3 columns named “start”, “end” and “globalState” encodes
the structural information. Each row provides the information of a
different structure. The number of CpG positions is given by
$end - start +1$, and each row start value needs to be as the previous
row end value + 1. The globalState encodes island and non-island
structures with the character value “U” for islands and “M” for
non-islands.

``` r
# Example 3 structures of length 100: non-island, island, non-island
infoStr <- data.frame(start = c(1, 101, 201),
                        end = c(100, 200, 300),
                        globalState = c("M", "U", "M"))
infoStr
#>   start end globalState
#> 1     1 100           M
#> 2   101 200           U
#> 3   201 300           M
```

``` r
# Example 2 islands: first of length 1, second of lenght 15
infoStr <- data.frame(start = c(1, 2),
                        end = c(1, 16),
                        globalState = c("U", "U"))
infoStr
#>   start end globalState
#> 1     1   1           U
#> 2     2  16           U
```

## Initial methylation state distribution

Can be sampled as specified in the section “Model parameters” when
infoStr has only start, end and globalState columns, or given by the
user as additional columns of the data frame. Note that when given, each
of the structures frequencies needs to sum up to 1.

``` r
# Example 3 structures of length 100 with customized initial methylation state equilibrium frequencies
infoStr <- data.frame(start = c(1, 101, 201),
                        end = c(100, 200, 300),
                        globalState = c("M", "U", "M"),
                        u_eqFreq = c(0.1, 0.8, 0.1),
                        p_eqFreq = c(0.1, 0.1, 0.1),
                        m_eqFreq = c(0.8, 0.1, 0.8))
infoStr
#>   start end globalState u_eqFreq p_eqFreq m_eqFreq
#> 1     1 100           M      0.1      0.1      0.8
#> 2   101 200           U      0.8      0.1      0.1
#> 3   201 300           M      0.1      0.1      0.8
```

## Parameter values

The default parameter values used by MethEvolRSIM can be obtained using
the function `get_model_parameters()`, that returns a data frame with
the columns named as specified in the section “Model parameters”:

``` r
default_paramValues <- get_parameterValues()
default_paramValues
#>   alpha_pI beta_pI alpha_mI beta_mI alpha_pNI beta_pNI alpha_mNI beta_mNI  mu
#> 1      0.1       1      0.1     0.5       0.1        1       0.5      0.1 0.1
#>   alpha_Ri iota
#> 1      0.1  0.3
```

The section “Model parameters” provides further information about each
of the parameters, including each parameter range of possible values. To
simulate data with a different parameterization, the new value of a
chosen parameter can be modified in the data frame.

``` r
# Example: modification of parameter iota to value 0.2
default_paramValues$iota <- 0.2
default_paramValues
#>   alpha_pI beta_pI alpha_mI beta_mI alpha_pNI beta_pNI alpha_mNI beta_mNI  mu
#> 1      0.1       1      0.1     0.5       0.1        1       0.5      0.1 0.1
#>   alpha_Ri iota
#> 1      0.1  0.2
```

## Initial DNA Methylation Data

### Simulation

The initial data is generated in an object of class
combiStructureGenerator.

``` r
# Example with customized initial methylation frequencies and customized parameter values
custom_infoStr <- data.frame(start = c(1, 101, 201),
                        end = c(100, 200, 300),
                        globalState = c("M", "U", "M"),
                        u_eqFreq = c(0.1, 0.8, 0.1),
                        p_eqFreq = c(0.1, 0.1, 0.1),
                        m_eqFreq = c(0.8, 0.1, 0.8))
custom_params <- get_parameterValues()
custom_params$mu <- 0.005
initialD <- simulate_initialData(infoStr = custom_infoStr, params = custom_params)
# Returns customized parameters and simulated data
initialD$params
#>   alpha_pI beta_pI alpha_mI beta_mI alpha_pNI beta_pNI alpha_mNI beta_mNI    mu
#> 1      0.1       1      0.1     0.5       0.1        1       0.5      0.1 0.005
#>   alpha_Ri iota
#> 1      0.1  0.3
initialD$data
#> <combiStructureGenerator>
#>   Public:
#>     add_offspring_index: function (i) 
#>     branch_evol: function (branch_length, dt, testing = FALSE) 
#>     clone: function (deep = FALSE) 
#>     copy: function () 
#>     get_island_index: function () 
#>     get_island_number: function () 
#>     get_IWE_events: function () 
#>     get_mu: function () 
#>     get_name: function () 
#>     get_offspring_index: function () 
#>     get_own_index: function () 
#>     get_parent_index: function () 
#>     get_singleStr: function (i) 
#>     get_singleStr_number: function () 
#>     initialize: function (infoStr, params = NULL, testing = FALSE) 
#>     set_IWE_events: function (a) 
#>     set_name: function (a) 
#>     set_offspring_index: function (i) 
#>     set_own_index: function (i) 
#>     set_parent_index: function (i) 
#>     set_singleStr: function (singStrList) 
#>   Private:
#>     interval_evol: function (interval_length, dt, testing = FALSE) 
#>     IWE_events: NULL
#>     IWE_rate: 0.005
#>     mu: 0.005
#>     name: NULL
#>     offspring_index: NULL
#>     own_index: NULL
#>     parent_index: NULL
#>     set_IWE_rate: function () 
#>     singleStr: list
#>     singleStr_end: 100 200 300
#>     singleStr_globalState: M U M
#>     singleStr_start: 1 101 201
#>     SSE_evol: function (dt, testing = FALSE)
```

``` r
# Example with simulated initial methylation frequencies and default parameter values
custom_infoStr <- data.frame(start = c(1, 101, 201),
                        end = c(100, 200, 300),
                        globalState = c("M", "U", "M"))
initialD <- simulate_initialData(infoStr = custom_infoStr)
# Returns default parameters
initialD$params
#>   alpha_pI beta_pI alpha_mI beta_mI alpha_pNI beta_pNI alpha_mNI beta_mNI  mu
#> 1      0.1       1      0.1     0.5       0.1        1       0.5      0.1 0.1
#>   alpha_Ri iota
#> 1      0.1  0.3
initialD$data
#> <combiStructureGenerator>
#>   Public:
#>     add_offspring_index: function (i) 
#>     branch_evol: function (branch_length, dt, testing = FALSE) 
#>     clone: function (deep = FALSE) 
#>     copy: function () 
#>     get_island_index: function () 
#>     get_island_number: function () 
#>     get_IWE_events: function () 
#>     get_mu: function () 
#>     get_name: function () 
#>     get_offspring_index: function () 
#>     get_own_index: function () 
#>     get_parent_index: function () 
#>     get_singleStr: function (i) 
#>     get_singleStr_number: function () 
#>     initialize: function (infoStr, params = NULL, testing = FALSE) 
#>     set_IWE_events: function (a) 
#>     set_name: function (a) 
#>     set_offspring_index: function (i) 
#>     set_own_index: function (i) 
#>     set_parent_index: function (i) 
#>     set_singleStr: function (singStrList) 
#>   Private:
#>     interval_evol: function (interval_length, dt, testing = FALSE) 
#>     IWE_events: NULL
#>     IWE_rate: 0.1
#>     mu: 0.1
#>     name: NULL
#>     offspring_index: NULL
#>     own_index: NULL
#>     parent_index: NULL
#>     set_IWE_rate: function () 
#>     singleStr: list
#>     singleStr_end: 100 200 300
#>     singleStr_globalState: M U M
#>     singleStr_start: 1 101 201
#>     SSE_evol: function (dt, testing = FALSE)
```

To get the parameter values of a combiStructureGenerator object, it can
be given as argument to the function `get_parameterValues()`

``` r
combiStr_object <- initialD$data
get_parameterValues(combiStr_object)
#>   alpha_pI beta_pI alpha_mI beta_mI alpha_pNI beta_pNI alpha_mNI beta_mNI  mu
#> 1      0.1       1      0.1     0.5       0.1        1       0.5      0.1 0.1
#>   alpha_Ri iota
#> 1      0.1  0.3
```

### Exploration of initial data

The methods of the class combiStructureGenerator can be used to access
each of the contained structures:

``` r
# E.g. access fist contained structure
combiStr_object$get_singleStr(1)
#> <singleStructureGenerator>
#>   Public:
#>     clone: function (deep = FALSE) 
#>     get_alpha_mI: function () 
#>     get_alpha_mNI: function () 
#>     get_alpha_pI: function () 
#>     get_alpha_pNI: function () 
#>     get_alpha_Ri: function () 
#>     get_beta_mI: function () 
#>     get_beta_mNI: function () 
#>     get_beta_pI: function () 
#>     get_beta_pNI: function () 
#>     get_combiStructure_index: function () 
#>     get_eqFreqs: function () 
#>     get_iota: function () 
#>     get_seq: function () 
#>     get_seq2ndButLastPos: function () 
#>     get_seq2ndPos: function () 
#>     get_seqFirstPos: function () 
#>     get_seqLastPos: function () 
#>     init_neighbSt: function () 
#>     initialize: function (globalState, n, eqFreqs = NULL, combiStr = NULL, combiStr_index = NULL, 
#>     initialize_ratetree: function () 
#>     IWE_evol: function (testing = FALSE) 
#>     private: environment
#>     self: singleStructureGenerator, R6
#>     SSE_evol: function (dt, testing = FALSE) 
#>     update_interStr_firstNeighbSt: function (leftNeighbSt, rightNeighbSt) 
#>     update_interStr_lastNeighbSt: function (leftNeighbSt, rightNeighbSt) 
#>   Private:
#>     alpha_mI: 0.1
#>     alpha_mNI: 0.5
#>     alpha_pI: 0.1
#>     alpha_pNI: 0.1
#>     alpha_Ri: 0.1
#>     beta_mI: 0.5
#>     beta_mNI: 0.1
#>     beta_pI: 1
#>     beta_pNI: 1
#>     check_ratetree: function (s) 
#>     choose_number_of_changes: function (dt) 
#>     choose_random_seqpos: function (testing = FALSE) 
#>     combiStructure_index: 1
#>     eqFreqs: 0.974429060485402 0.0090046575021602 0.0165662820124377
#>     get_leftStr_neighbSt: function () 
#>     get_nextStr: function () 
#>     get_prevStr: function () 
#>     get_rightStr_neighbSt: function () 
#>     globalState: M
#>     init_Rc_values: function () 
#>     init_Ri_values: function () 
#>     iota: 0.3
#>     mapNeighbSt_matrix: 1 4 7 2 5 8 3 6 9
#>     my_combiStructure: combiStructureGenerator, R6
#>     neighbSt: 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1  ...
#>     Q: list
#>     Qc: list
#>     Qi: list
#>     ratetree: list
#>     Rc_values: list
#>     Ri_values: 2.80505952795888e-06 0.00577075223750925 0.894226446338048
#>     sample_eqFreqs: function () 
#>     seq: 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1  ...
#>     set_Q: function () 
#>     set_Qc: function () 
#>     set_Qi: function () 
#>     siteR: 3 2 2 2 1 3 2 1 1 1 2 3 3 3 3 3 1 2 2 1 1 3 2 2 1 2 1 3  ...
#>     update_intraStr_neighbSt: function (position) 
#>     update_neighbSt: function (position) 
#>     update_ratetree: function (position, rate)
```

Each contained structure is an object of singleStructureGenerator and
its methods can be used to get the contained information.

``` r
# E.g. get the equilibrium frequencies of the first singleStructureGenerator
combiStr_object$get_singleStr(1)$get_eqFreqs()
#> [1] 0.974429060 0.009004658 0.016566282
```

``` r
# E.g. get the sequence of methylation states of the second singleStructureGenerator
combiStr_object$get_singleStr(1)$get_seq()
#>   [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#>  [38] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 1
#>  [75] 1 1 1 1 1 1 3 1 1 1 1 1 1 1 1 1 1 1 1 1 1 3 1 1 1 2
```

## Simulate DNA Methylation Data Evolution

The evolution data is generated in an object of class
treeMultiRegionSimulator.

``` r
# Example with customized initial methylation frequencies, customized parameter values and default dt (0.01)
tree <- "((a:1, b:1):2, c:2, (d:3.7, (e:4, f:1):3):5);"
custom_infoStr <- data.frame(start = c(1, 101, 201),
                        end = c(100, 200, 300),
                        globalState = c("M", "U", "M"),
                        u_eqFreq = c(0.1, 0.8, 0.1),
                        p_eqFreq = c(0.1, 0.1, 0.1),
                        m_eqFreq = c(0.8, 0.1, 0.8))
custom_params <- get_parameterValues()
custom_params$mu <- 0.005
evolD <- simulate_evolData(infoStr = custom_infoStr, tree = tree, params = custom_params)
#> [1] "Simulating data at root and letting it evolve along given tree:  ((a:1, b:1):2, c:2, (d:3.7, (e:4, f:1):3):5);"
#> All validation states are valid. No errors found.
#> All validation states are valid. No errors found.
#> All validation states are valid. No errors found.
#> All validation states are valid. No errors found.
#> All validation states are valid. No errors found.
#> All validation states are valid. No errors found.
#> All validation states are valid. No errors found.
# Returns customized parameters, tree used, time step length for SSE process used (dt) and simulated data
evolD$data
#> <treeMultiRegionSimulator>
#>   Public:
#>     Branch: list
#>     branchLength: NA 2 1 1 2 5 3.7 3 4 1
#>     clone: function (deep = FALSE) 
#>     initialize: function (infoStr = NULL, rootData = NULL, tree, params = NULL, 
#>     treeEvol: function (T, dt = 0.01, parent_index = 1, testing = FALSE)
evolD$dt
#> [1] 0.01
evolD$tree
#> [1] "((a:1, b:1):2, c:2, (d:3.7, (e:4, f:1):3):5);"
evolD$params
#>   alpha_pI beta_pI alpha_mI beta_mI alpha_pNI beta_pNI alpha_mNI beta_mNI    mu
#> 1      0.1       1      0.1     0.5       0.1        1       0.5      0.1 0.005
#>   alpha_Ri iota
#> 1      0.1  0.3
```

If the evolution of given data at the root is simulated, and it wants to
be done with customized parameter values, then the customized parameters
are given to `simulate_initialData()` instead to `simulate_evolData()`.

``` r
# Example with customized initial methylation frequencies, customized parameter values and default dt (0.01)
tree <- "((a:1, b:1):2, c:2, (d:3.7, (e:4, f:1):3):5);"
custom_infoStr <- data.frame(start = c(1, 101, 201),
                        end = c(100, 200, 300),
                        globalState = c("M", "U", "M"),
                        u_eqFreq = c(0.1, 0.8, 0.1),
                        p_eqFreq = c(0.1, 0.1, 0.1),
                        m_eqFreq = c(0.8, 0.1, 0.8))
custom_params <- get_parameterValues()
custom_params$mu <- 0.005
initialD <- simulate_initialData(infoStr = custom_infoStr, params = custom_params)$data
evolD <- simulate_evolData(rootData =initialD, tree = tree)
#> [1] "Parameter values set as in given rootData"
#> [1] "Simulating evolution of given data at root along given tree:  ((a:1, b:1):2, c:2, (d:3.7, (e:4, f:1):3):5);"
#> All validation states are valid. No errors found.
#> All validation states are valid. No errors found.
#> All validation states are valid. No errors found.
#> All validation states are valid. No errors found.
#> All validation states are valid. No errors found.
#> All validation states are valid. No errors found.
#> All validation states are valid. No errors found.
```

### Exploration of evolved data

The treeMultiStructureGenerator object contains a list `$Branch` with
the combiStructureGenerator representing the data at the each node after
a tree branch (but for the data at the root).

To know which offspring branches a node has, the combiStructureGenerator
object method `$get_offspring_index()` can be used.

``` r
treeData <- evolD$data
# E.g. offspring of data at root 
treeData$Branch[[1]]$get_offspring_index()
#> [1] 2 5 6
```

If the node is an inner node, its name is set as NULL. It can be
accessed with the combiStructureGenerator method `$get_name()`

``` r
treeData$Branch[[1]]$get_name()
#> NULL
```

But the leafs are named as in the given tree.

``` r
treeData$Branch[[4]]$get_name()
#> [1] "b"
```

If there were IWE events sampled in a branch, it can be accessed with
the combiStructureGenerator method `$get_IWE_events()`

``` r
treeData$Branch[[4]]$get_IWE_events()
#> [1] FALSE
```

In the case we want to know the equilibrium frequencies of a
singleStructureGenerator object contained in the combiStructureGenerator
we can use the combiStructureGenerator method to access the
singleStructureGenerator object of interest `get_singleStr(i)` followed
by the singleStructureGenerator method `get_eqFreqs()`

``` r
# E.g. get the ewuilibrium frequencies of the second structure at the leaf "b"
treeData$Branch[[4]]$get_singleStr(2)$get_eqFreqs()
#> [1] 0.8 0.1 0.1
```

# Model parameters

## Initial methylation state distribution

Structures are categorized as CpG islands or non-islands. In the root of
the genealogy we assume that each structure samples an equilibrium
probability triple $(\pi_u,\pi_p, \pi_m)$. Current parameter values
favor prevalent homogeneity in methylation levels within structures,
most CpG islands with a majority predominantly unmethylated and a small
proportion predominantly methylated, and most non-islands with a
majority predominantly methylated and a small proportion predominantly
unmethylated. The sampling of the equilibrium probability tripple is
done using two distinct Beta distributions to sample the corresponding
methylation frequencies.

**Island parameterization: alpha_pI, beta_pI, alpha_mI, beta_mI**

MethEvolRSIM samples $\pi_p$ from the first Beta distribution. Current
parameterization uses as initial parameter values a relatively small
$\alpha_{I_p}$ and relatively large $\beta_{I_p}$ to encourage minor
assignment of partially methylated states.

The second Beta distribution, scaled by $1 - \pi_p$, corresponds to the
proportion of methylated and unmethylated states. Currently both
parameter values ($\alpha_{I_m}$ and $\beta_{I_m}$) set to small values
while ${\alpha_{I_m} < \beta_{I_m}}$ to favour shifting the stochastic
choice toward homogeneous states with a higher proportion of
predominantly unmethylated islands, and assigned to $\pi_m$. Finally, we
set $\pi_u$ as $1 - \pi_p - \pi_m$.

**Non-island parameterization: alpha_pNI, beta_pNI, alpha_mNI,
beta_mNI**

Following the same rationale as for island parameterization, all the
regions not belonging to an island sample a common equilibrium
probability triple from two Beta distributions with parameters
$\alpha_{NI_p}$, $\beta_{NI_p}$, $\alpha_{NI_m}$ and $\beta_{NI_m}$. The
parameter values are currently set to favor the characteristic high
overall methylation level.

You can simulate DNA methylation along a tree given some data at root,
or by simulating data at root given a structural distribution and
letting it evolve.

**Parameter ranges** The alpha and beta parameters of a Beta
distribution must be greater than 0, that is (0, $\infty$).

## DNA methylation evolution processes

DNA methylation state evolution happens by two structure-specific
events: CpG single site events (SSEs) and CpG island wide events (IWEs).

### Island-wide events (IWEs)

IWEs change the island methylation probabilities and some of the CpG
sites in the CpG island can simultaneously change their state at the
same time. Different CpG islands, however, are assumed to evolve
independently of each other. IWEs occur per CpG island at a rate denoted
by $\mu$. Consequently, over a branch of length $l$, the expected number
of IWEs per island is given by $\mu \cdot l$.

**IWE parameterization: mu. Parameter range:** Should be a small number,
but it can take any value from 0, that is \[0, $\infty$)

### Single-site events (SSEs)

In addition to IWEs we allow for single site events (SSEs), which change
the methylation states of single CpG sites within and between CpG
islands. In our model SSEs account for two dynamics: independent of each
other (SSEi), and assuming correlations between adjacent sites (SSEc).

As SSEc transition probabilities are dependent on adjacent sites we need
to approximate SSE based evolution by short time steps, allowing for CpG
sites to actualize their state. Therefore, we first discretize the
branch intervals between IWEs in time intervals of length 0.01. For each
CpG position and each time interval the rate matrix for the transitions
between the states unmethylated, partially methylated and methylated
results from the addition of the rate matrices of the two types of SSE
events, so that $\mathbb{E}(R_i + R_c)=1$, $R_i$ and $R_c$ representing
the rate factors of independent and correlated processes. The proportion
of each of the two SSE types is determined by the parameter $\iota$.

Different positions $j$ have different SSEi rates $R_{i,j}$ coming from
a discretized gamma distribution with 3 categories, expected value
$\iota$ and shape parameters $\alpha_{R_{i}}$ and
${\beta_{R_{i}}}=\frac{\alpha_{R_{i}}}{\iota}$. The probability to be in
each respective rate category is 1/3. For the rate of collaborative SSEs
$R_c = 1-\iota$, the probability of considering left of right neighbour
is equal.

**SSE process parameterization: iota, alpha_Ri** **Parameter ranges**
The alpha parameter of a Gamma distribution must be greater than 0, that
is (0, $\infty$). Iota must be between 0 (non-included) and 1, that is
(0,1\].
