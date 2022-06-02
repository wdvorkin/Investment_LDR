# Multi-Stage Investment Decision Rules for Power Systems with Performance Guarantees

This repository contains supplementary materials, optimization data, and codes for [Multi-Stage Investment Decision Rules for Power Systems with Performance Guarantees](arxiv.org) by [Vladimir Dvorkin](http://wdvorkin.github.io), [Dharik Mallapragada](https://mallapragada.mit.edu), and [Audun Botterud](https://botterud.mit.edu).

<img width="1082" alt="Screen Shot 2022-06-02 at 9 53 04 AM" src="https://user-images.githubusercontent.com/31773955/171645236-63fd26b1-8419-4273-a55b-e8f957e14f7f.png">

The optimization dataset is an ensemble of M.Sc. [thesis](https://dspace.mit.edu/bitstream/handle/1721.1/140416/schwartz-aaronms-sm-tpp-2021.pdf?sequence=1&isAllowed=y) of Aaron Schwartz, the [annual technology baseline](https://atb.nrel.gov/electricity/2021/index) and [future electrification study](https://www.nrel.gov/docs/fy21osti/72330.pdf) reports by the National Renewable Energy Laboratory (NREL), and the [world energy outlook 2021](https://www.iea.org/reports/world-energy-outlook-2021) by the International Energy Agency. The authors acknowledge Jack Moris's effort in assembling the data. 

Please, refer to [Apendix.pdf](https://github.com/wdvorkin/Investment_LDR/blob/main/Appendix.pdf) for supplementary materials. 

The optimization models were implemented in Julia language (v.1.6) using JuMP modeling language for mathematical optimization embedded in Julia. The models require Mosek comercial optimization solver, which needs to be installed and licensed.

To activate the packages in ```Project.toml```, clone the project using e.g. ```git clone```, ```cd``` to the project directory and call
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

<!-- By default, the program returns the solution of the ```"CC_OPF"``` mechanism and stores the results in ```~/output/CC_OPF```. To run the other mechanisms, parse ```"D_OPF"```, ```"ToV_CC_OPF"```, ```"TaV_CC_OPF"``` or ```"CVaR_CC_OPF"``` using option ```-m```, e.g. 
```
$ julia DP_CC_OPF.jl -m "CVaR_CC_OPF"
```
The results will be stored in ```~/output/CVaR_CC_OPF```.  -->


