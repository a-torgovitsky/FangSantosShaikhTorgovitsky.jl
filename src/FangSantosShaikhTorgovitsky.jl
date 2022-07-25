module FangSantosShaikhTorgovitsky
    using CSV
    using CategoricalArrays
    using DataFrames
    using DelimitedFiles
    using Distributed
    using Distributions
    using GLM
    using Gurobi
    using JuMP
    using LinearAlgebra
    using MathOptFormat
    using MathOptInterface
    using Printf
    using Random
    using Roots
    using Parameters
    using Sobol
    using Statistics
    using StatsBase
    using PrettyTables
    using LoggingExtras
    using Dates

    include("menu.jl")
    include("lpmodels.jl")
    include("cone.jl")
    include("gmm.jl")
    include("mixedlogit.jl")
    include("montecarlo.jl")
    include("default_globals.jl")
    include("rundesigns.jl")
    include("paper_simulations.jl")
    include("util_opt.jl")
    include("util_print.jl")
    include("util_testing.jl")
    include("util_save.jl")

    # https://github.com/jump-dev/Gurobi.jl/issues/388
    # Define one Gurobi environment and keep reusing it
    const GLOBAL_GENV = Ref{Gurobi.Env}()
    function __init__()
        global GLOBAL_GENV[] = Gurobi.Env()
    end
end
