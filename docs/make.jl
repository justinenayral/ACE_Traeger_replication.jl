
push!(LOAD_PATH,"../src/")


using ACE_Traeger_replication

using Documenter


#using DocumenterTools
makedocs(
         sitename = "ACE_Traeger_replication.jl",
         modules  = [ACE_Traeger_replication],
         pages=[
                "Home" => "index.md"
               ],
               format = Documenter.HTML(prettyurls = false)
)


deploydocs(
    repo="github.com:justinenayral/ACE_Traeger_replication.jl.git",
    devbranch="main"
)


