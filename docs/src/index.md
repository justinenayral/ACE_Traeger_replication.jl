# ACE_Traeger_Replication.jl

**Authors**: Norbert Monti, Nayral Justine


This package reproduces the findings of Traeger, Christian P. (2023) in his paper titled 'ACE—Analytic Climate Economy,' published in the *American Economic Journal: Economic Policy*, Volume 15, Issue 3, pages 372-406. While the [original replication materials] (https://www.openicpsr.org/openicpsr/project/154141/version/V1/view) provided by the authors were coded in Matlab, we have used Julia to build a replication package, encompassing the main results, Figure II, Figure III, Figure IV, and Table I.

The paper examines optimal carbon taxation using integrated assessment models (IAMs) of climate change. These models are designed to evaluate the long-term interactions among economic production, greenhouse gas emissions, and global warming. C. Traeger discusses the implications of temperature and carbon tax impact. The persistence of carbon increases the optimal tax twofold to thirtyfold, depending on the calibration. On the contrary, the delay in temperature dynamics (Ocean cooling) decreases the carbon tax from 65 to 25%. The Analytics Climate Economy (ACE) model is closed to Nordhaus' DICE model. It incorporate most elements of IAMs. Labor, capital, technology and energy produced output are either consumed or invested. The author distinguish "Dirty" energy sectors, consuming fossil fuels and generating greenhouse gas emissions. These gases accumulate in the atmorsphere causing radiative forcing and increase global temperature, which reduces output. This economic model aims at helping economists to develop more accurate opinions about the social cost of carbon.

## Data availability:
Our replication packages require to download the data used by the authors from the orginal [replication package](https://www.openicpsr.org/openicpsr/project/154141/version/V1/view?flag=follow&pageSize=100&sortOrder=(?title)&sortAsc=true).

## Path:
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
