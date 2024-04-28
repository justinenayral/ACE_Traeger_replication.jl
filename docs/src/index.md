# ACE_Traeger_Replication.jl

**Authors**: Norbert Monti, Nayral Justine


This package reproduces the findings of Traeger, Christian P. (2023) in his paper titled 'ACE—Analytic Climate Economy,' published in the *American Economic Journal: Economic Policy*, Volume 15, Issue 3, pages 372-406. While the [original replication materials] (https://www.openicpsr.org/openicpsr/project/154141/version/V1/view) provided by the authors were coded in Matlab, we have used Julia to build a replication package, encompassing the main results, Figure II, Figure III, Figure IV, and Table I.

The paper examines optimal carbon taxation using integrated assessment models (IAMs) of climate change. These models are designed to evaluate the long-term interactions among economic production, greenhouse gas emissions, and global warming. C. Traeger discusses the implications of temperature and carbon tax impact. The persistence of carbon increases the optimal tax twofold to thirtyfold, depending on the calibration. On the contrary, the delay in temperature dynamics (Ocean cooling) decreases the carbon tax from 65 to 25%. The Analytics Climate Economy (ACE) model is closed to Nordhaus' DICE model. It incorporate most elements of IAMs. Labor, capital, technology and energy produced output are either consumed or invested. The author distinguish "Dirty" energy sectors, consuming fossil fuels and generating greenhouse gas emissions. These gases accumulate in the atmorsphere causing radiative forcing and increase global temperature, which reduces output. This economic model aims at helping economists to develop more accurate opinions about the social cost of carbon.

## Data availability
Our replication packages require to download the data used by the authors from the orginal [replication package](https://www.openicpsr.org/openicpsr/project/154141/version/V1/view?flag=follow&pageSize=100&sortOrder=(?title)&sortAsc=true).

## Path
The user needs to save the datas and the julia package in path:
```julia
path = "C:/Users/..."
```

Figures generated are saved in graph path:
```julia
graphpath=path*"figure/"
```
## Packages
Before using the package, it is required to install:
```julia
using Pkg 
Pkg.add("MAT") #To open .mat datasets
Pkg.add("Printf") #For labelling graphs
Pkg.add("LinearAlgebra")
Pkg.add("DataFrames")
Pkg.add("CSV")
Pkg.add("XLSX")
Pkg.add("NLsolve")
Pkg.add(("Plots")
```
## Dammage functions
The authors used differents definition for calculating environment damages. This section presents the diferent functions the user can use to calculate them, depending of the method he prefers. In case of issues, while running the code just use the ? in the repl, followed by the function's name, to get the documentation.

#### ACE
The damage function in ACE takes an exponential functional form: 
``D(T_{1,t})=1-\exp(-\xi_0\exp(\xi_1 T_{1,t})+\xi_0)``
with ``\xi_0``, a free damage parameters and the paramteter ``\xi_1=\frac{\log 2}{s}`` pinned down by climate sensitivity *s*.

*Input:*
- ``\xi_0``: xi0 (Calibrated to 0.022 in the paper. The base calibration is an exact match of the two calibration points 0° and 2.5° in the 2007 model)
- ``\xi_1``: xi1 (Estimated to ``\frac{1}{4}`` in the paper)
- ``T``: temperature

*Syntax:*
```julia
dam_ACE(xi0, xi1, T)
```
*Output:*
The calculated ACE damage estimate.

#### DICE
DICE damage function replace the quadratic damage term ``aT^2``, to limit damages to 100 percent of production. 
``D(T)=1-\frac{1}{1+aT^2}``

*Input:*
- ``a`` : Damage coefficients (``a=0.0028`` for 2007, ``a=0.00267`` for 2013, ``a=0.00236`` for 2016)
- ``T``: temperature

*Syntax:*
```julia
ACE_Traeger_replication.dam_DICE(a, T)
```
*Output:*
The calculated DICE damage estimate.

#### Howard and Sterner (2017)
This damage function is based on Howard and Sterner (2017):
``D(T)=aT^2``

*Input:*
- ``a`` : Damage coefficient. By default, it is set to 0.01145. The user must write the value of a explicitely, if he wants to use the dam_Sterner function with another ``a``.
- ``T``: temperature

*Syntax:*
```julia
ACE_Traeger_replication.dam_Sterner(T, a)
```

*Output:*
The calculated DICE-Howard-Sterner damage estimate.

## Figure II
The function below replicated original figure II in the paper. This figure shows the predicted damages of the different model depending on the temperature degrees above the preindustrial level.
The different models are:
- DICE 2007, 2013, 2016
- ACE base damage calibration. It matches the two calibration points 0 and 2.5° in the 2007 model. 
- HSP-norm: Howard and Sterner's damage function using DICE's approach to limit damages to 100 percents
- HSP-nn: Howard and Sterner's damage function using DICE's approach, not normalized to limit damages to 100 percents. Damages exceed production at 9.5°. 

The left side of the plot refers to the temperature range of the IPCC secenarios in Figure 3 and the right zone is just a focus on lower degrees of warming. 

*Input:*
- path

*Syntax:*
```julia
ACE_Traeger_replication.figure2(path)
```

*Output:*
- Figure 2 saved into graphpath

Original figure: 
![](figure2or.png)

The figure generated by our package: 
![](figure2.png)


## Figure III
*Input:*
- path

*Syntax:*
```julia
ACE_Traeger_replication.figure3(path)
```

*Output:*
- Figure 3 saved into graphpath

Original figure: 
![](figure3or.png)

The figure generated by our package: 
![](figure3.png)
