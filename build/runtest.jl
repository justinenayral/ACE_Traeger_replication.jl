module runtest

export @testset

using MAT
using Printf
using LinearAlgebra
using NLsolve
using DataFrames
using Plots
using CSV
using XLSX
using Test

include("ACE_Traeger_replication.jl")


@testset "Damage function ACE" begin
    @test abs(ACE_Traeger_replication.dam_ACE(0.022, 0.25, 3)-0.02427451779)<0.000000001
end

@testset "Damage function DICE" begin
    @test abs(ACE_Traeger_replication.dam_DICE(0.0028,3)-0.02458056964494726)<0.000000001
end

@testset "Damage function ACE" begin
    @test abs(ACE_Traeger_replication.dam_Sterner(3)-0.10305)<0.000000001
end

end