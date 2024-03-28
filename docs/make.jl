using OptTestFun
using Documenter

DocMeta.setdocmeta!(OptTestFun, :DocTestSetup, :(using OptTestFun); recursive=true)

makedocs(;
    modules=[OptTestFun],
    authors="AllanAmorim <150050564+AllanAmorim@users.noreply.github.com> and contributors",
    sitename="OptTestFun.jl",
    format=Documenter.HTML(;
        canonical="https://AllanAmorim.github.io/OptTestFun.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/AllanAmorim/OptTestFun.jl",
    devbranch="master",
)
