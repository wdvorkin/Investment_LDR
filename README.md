# Multi-Stage Investment Decision Rules for Power Systems with Performance Guarantees

This repository contains supplementary materials, optimization data, and codes for [Multi-Stage Investment Decision Rules for Power Systems with Performance Guarantees](arxiv.org) by [Vladimir Dvorkin](http://wdvorkin.github.io), [Dharik Mallapragada](https://mallapragada.mit.edu), and [Audun Botterud](https://botterud.mit.edu).

<img width="1082" alt="Screen Shot 2022-06-02 at 9 53 04 AM" src="https://user-images.githubusercontent.com/31773955/171645236-63fd26b1-8419-4273-a55b-e8f957e14f7f.png">

The optimization dataset is an ensemble of M.Sc. [thesis](https://dspace.mit.edu/bitstream/handle/1721.1/140416/schwartz-aaronms-sm-tpp-2021.pdf?sequence=1&isAllowed=y) of Aaron Schwartz, the [annual technology baseline](https://atb.nrel.gov/electricity/2021/index) and [future electrification study](https://www.nrel.gov/docs/fy21osti/72330.pdf) reports by the National Renewable Energy Laboratory (NREL), and the [world energy outlook 2021](https://www.iea.org/reports/world-energy-outlook-2021) by the International Energy Agency. The authors acknowledge Jack Moris's effort in assembling the data. 

Please, refer to [Apendix.pdf](https://github.com/wdvorkin/Investment_LDR/files/8832872/appendix.pdf) for supplementary materials. 

The optimization models were implemented in Julia language (v.1.6) using [JuMP](https://github.com/jump-dev/JuMP.jl) modeling language for mathematical optimization and commercial [Mosek](https://github.com/MOSEK/Mosek.jl) optimization solver, both embedded in Julia. The solver needs to be licensed (free for academic use).

To activate and run the project, clone this repository, e.g., using ```git clone```, then ```cd``` to the project directory and call
```
$ julia 
julia> ]
(@v1.6) pkg> activate .
(Investment_LDR) pkg> instantiate
```

where ```julia``` is an alias to Julia installation. To run the code, ```cd``` to the project directory and call
```
$ julia main.jl
```

By default, the code returns the solution of the multi-stage investment planning with chance constraints reformulation under the Normal distribution assumption. The results will be stored in the ```output``` folder in the project root. To solve the planning problem in a distributionally robust manner, set option ```-r``` to ```DRO-DS```, where ```DS``` stands for double-sided chance constraint reformulation, i.e., 
```
$ julia main.jl -r "DRO-DS"
```
To see all available options, type 
```
$ julia main.jl --help
```

