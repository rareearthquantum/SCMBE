
# Getting Started

## Step 1: Prerequisites

Download Julia and Git.

## Step 2: Activating the Project and Downloading Dependencies

```
cd whereyouwantit
git clone https://github.com/rareearthquantum/SCMBE.git
cd SCMBE
julia --project=.
```
Then within the just opened Julia REPL, enter Pkg mode by typing
```
]
```
Then download all the necessary dependencies by typing
```
instantiate
```

Exit the Julia REPL using `CTRL+D`.

## Step 3: Running an Example Script

Following directly from the previous step, to run an example script enter in
```
julia --project=. examples/IMPORTANT_EXAMPLE.jl
```