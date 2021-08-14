# Tensor-Networks

## 2D Ising Model

### To do:

- [ ] Define a function to generate Ising Model:
  - [ ] initialise network
  - [ ] calculate partition function (perform partial trace on each site $\text{Tr}\left(\text{e}^{-\beta H}\right)$ (T -> O))
- [ ] Write a CTMRG (Corner Transfer Matrix Renormalisation Group) algorithm:
  - [ ] define corner tensors 
  - [ ] define edge tensors
- [ ] Perform CTRMG:
  - [ ] expectation value: $\frac{\text{Tr}\left(\text{e}^{-\beta H} \sigma_z\right)}{\text{Tr}\left(\text{e}^{-\beta H}\right)}$

### Questions:

1. What are the dimesnions of the network?
2. What is the dimenion of chi?
3. How te replace infinite parts of the network with finite elements?
4. Are boundary conditions in the model finite or infinite?
5. Are all the corner or edge tensors the same?
6. When performing CTMRG with sigma_z, should sigma_z be considered in edge and corner tensors?