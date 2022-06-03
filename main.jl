#!/usr/bin/env julia
# using Pkg
# Pkg.activate("../.")
# Pkg.instantiate()

using Statistics, LinearAlgebra, Distributions
using JuMP, Mosek, MosekTools
using DataFrames, ArgParse, FileIO, JSON, CSV

# parse arguments
cd(dirname(@__FILE__))
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--emission", "-e"
            help = "probability of emission cap constraint violation probability"
            arg_type = Float64
            default = 0.1
        "--conviol", "-c"
            help = "probability of individual constraint violation"
            arg_type = Float64
            default = 0.001
        "--varcoef", "-v"
            help = "investment decision variance constraint"
            arg_type = Float64
            default = 10.0
        "--sigma", "-s"
            help = "standard deviation σ of ξ"
            arg_type = Float64
            default = 0.25
        "--samsize", "-n"
            help = "sample size of out-of-sample analysis"
            arg_type = Int64
            default = 1000
        "--dist", "-d"
            help = "distribution for out-of-sample analysis: Normal, Uniform, Laplace, Logistic"
            arg_type = String
            default = "Normal"
        "--ref", "-r"
            help = "reformulation of chance constraints: Normal or DRO-SS or DRO-DS"
            arg_type = String
            default = "Normal"
    end
    return parse_args(s)
end
args = parse_commandline()
outdir = "output_var_$(args["varcoef"])_sigma_$(args["sigma"])_emiss_$(args["emission"])_ref_$(args["ref"])_dist_$(args["dist"])/"
mkpath(outdir)

# load functions
include("scr/functions.jl")

# Load network data
caseID="data/southeast_US"
# load data
set = Dict(
:T => 5, :n_t => [1 4 7 10 13], :H => 24, :W => 14,
# load rate of change
:r_l => [1 0.08 0.175 0.175 0.175],
# capex rate of change
:r_q_wind => [1	-0.189	-0.041	-0.041	-0.041],
:r_q_pv   => [1	-0.279	-0.032	-0.032	-0.032],
:r_q_nucl => [1	-0.028	-0.034	-0.033	-0.031],
:r_q_ng   => [1	-0.026	-0.021	-0.021	-0.019],
:r_q_s_e   => [1	-0.224	-0.048	-0.048	-0.048],
:r_q_s_p   => [1	-0.105	-0.056	-0.056	-0.056],
# oper and maintance fixed cost rate of change
:r_o_wind => [1	-0.088	-0.059	-0.044	-0.036],
:r_o_pv   => [1	-0.168	-0.021	-0.021	-0.021],
:r_o_nucl => [1	-0.000	-0.000	-0.000	-0.000],
:r_o_ng   => [1	-0.000	-0.000	-0.000	-0.000],
:r_o_s    => [1	-0.000	-0.000	-0.000	-0.000],
# fuel cost rate of change
:r_ng => [1	0.11	0.06	0.01	-0.01],
:r_cl => [1	-0.033	-0.012	0.000	-0.009],
:year_per_period => 5, :int_rate => 0.045,
:cost => "linear", :var_pen => args["varcoef"],
:ε̅ => 0.01, :ε̅_e => args["emission"], :ε̅_f => 0.15, :ε̅_g => 0.05, :ε̅_s => 0.05, :ε̅_i => 0.1,
:α_y => args["varcoef"], :α_ϑ =>  args["varcoef"], :α_φ =>  args["varcoef"],
:β => NaN,
# investment technology limits
:y̅_wd => [79000 79000 79000 1000000 1000000],
:y̅_pv => [48000 48000 48000 1000000 1000000],
:e̅ => [150 125 100 75 50],
# uncertainty data
:σ => args["sigma"],
:dist => args["dist"],
# uncertainty modeling
:Ns => args["samsize"],
:cc_ref => args["ref"],
:cc_ref_double_sided => false,
:VoLL => true
)

# load data
data = load_data(caseID,set)
# load forecast
ξ = error_stat(set)

# solve primal models
sol_det_prim = determenistic_expansion_prim(set,data,ξ)
sol_sto_prim = stochastic_expansion_prim(data,set,ξ)

# # uncomment solve dual models
# sol_det_dual = determenistic_expansion_dual(set,data,ξ)
# sol_sto_dual = stochastic_expansion_dual(data,set,ξ)

# run out-of-sample analysis and estiamate constraint violation
post_processing(outdir,data,set,ξ)
