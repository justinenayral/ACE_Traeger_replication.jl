push!(LOAD_PATH,"../src/")
#using ACE_Traeger_replication
using Documenter
makedocs(
         sitename = "ACE_Traeger_replication.jl",
         modules  = [ACE_Traeger_replication],
         pages=[
                "Home" => "index.md"
               ]
)

deploydocs(
    repo="github.com:justinenayral/ACE_Traeger_replication.jl.git"
)