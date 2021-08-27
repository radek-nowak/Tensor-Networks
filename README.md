# Tensor-Networks

### To do:

## 2D Ising Model

### To do:

- [x] Define a function to generate Ising Model:
  - [x] initialise network
  - [x] calculate partition function (perform partial trace on each site)
 (T -> O))
- [ ] Write a CTMRG (Corner Transfer Matrix Renormalisation Group) algorithm:
  - [x] define corner tensors 
  - [x] define edge tensors
- [ ] Perform CTRMG:
  - [ ] calculate expectation value of `sigma_z`

### Questions:

1. What are the dimesnions of the network?
2. What is the dimenion of chi?
3. How to replace infinite parts of the network with finite elements (is it just an approximation?)?
4. Are boundary conditions in the model finite or infinite?
5. Are all the corner or edge tensors the same?
6. When performing CTMRG with `sigma_z`, should `sigma_z` be considered in edge and corner tensors?
7. How to perform partial trace on net T?


