#1 MVP
* create AST
* represent and print the example program using that AST

#2 Optimizations
* create an inspector to traverse input matrix and generate reachset
* create an AST traversal to prune iteration space based on reachset (VIPrune) 

# Thursday Jan 14 
## need to setup test for  making sure the solver is correct
Running the test, I run into an issue with the sparse solve. Try a different one?
The other solver won't work on sparse systems. Trying the other matrix,
also failed.
Try converting to array first: that works, but later it complains
about A not being triangular. I guess that's because it is, so I
need to slice off the upper part for the solver to be happy.
tril() complains about missing a positional argument.
Once again, shape is empty. Trying to toarray trick. 
Now it complains about the matrix being singular.

This is taking too much time. Instead, I will verify the solution
by running the c++ solver, and using scipy to multiply Ax and verify
that it is equal to b.

# Friday, Jan 15
This is taking too much time. I need to prioritize. 
For now, a morphism is considered "correct" iff it
produces the same output as the naive implementation.

# Saturday, Jan 16
Got pruning kind-of working.
Important: reachset needs to be ordered!!!!!!!!


