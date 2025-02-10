# A version of the standard GTAP model (V7.1) coded in Julia

The code folder contains three files:
1. GTAP0.jl--an initial attempt that uses numeric indices.
2. GTAP.ipynb--A Jupyter file with the model split into key components. IT is easy to add simulations at the end of the file and will eventually link to Julia's plotting ability.
3. GTAP.jl--Virtually the same code as in GTAP.ipynb.

The model comes with two sample databases. Both are sourced from GTAP V9 RC2 and are freely available. The first, 3x3, is a very small database with no land or natural resources. The second, is called DBUG0, and has somewhat more structure--essentially a 10x10 database, but also using GTAP's power database. The make matrix is non-diagonal and the different electricity activities collapse to one commodity. Users can choose which database at the beginning of either file.

At this stage, the model is set to mostly replicate V7.1 of the GTAP model coded in GEMPACK. It lacks a few features, such as an alternative closure for the allocation of global savings. (It does have two closures not in the GTAP model--fixed capital flows, or fixed capital flows as a share of GDP. The latter requires a residual region.) Next iterations will include energy and emissions. The model code can also be simplified, for example collapsing all Armington agents into a single set. And eventually, we may wish to integrate some form of flexible nesting. Finally, we will also want to link it seemlessly to an aggregation facility of the GTAP database.
