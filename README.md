
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
then within the Julia REPL enter in
```
]
```
which enters the Pkg mode, then enter
```
instantiate
```
to download all the necessary dependencies.

To exit the Julia REPL, use `CTRL+D`.

## Step 3: Running an Example Script

Following directly from the previous step, to run an example script enter in
```
julia --project=. examples/IMPORTANT_EXAMPLE.jl
```