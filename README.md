# FWI
Acoustic Full-waveform inversion (FWI) for isotropic media, where the wave equation is solved in the frequency domain with perfectly mathced layers (PML) and absorbing boundary condition (ABC).\
Forward and adjoint wavefields are cross-correlated to calculate the gradient for velocity correction.

# Generate true recordings
Run 'generate_true_rec.m' to generate true recordings 'true_rect.mat' based on a toy model 'trueModel'. Below is our true velocity model (stolen from the courst TPG4155 at NTNU).
![true velocity](https://user-images.githubusercontent.com/45905048/68906553-f81c5b80-0744-11ea-90ab-1384d10d7f28.jpg)
# Run FWI with true recordings
Run 'FWI.mat' to process FWI based on the true recordings to reconstruct the velocity model. The starting model is given below.
![starting velocity](https://user-images.githubusercontent.com/45905048/68906711-855fb000-0745-11ea-812b-d576b0eb66ec.jpg)
Below is the final result after 30 times of interations. The final velocity model is close to our true velocity model. The salt structure on the left of the true model is not well reconstructed. The reason is that the high velocity of the salt bends the P-wave considerably, making it hard for P-wave penetrate the inner side. 
![progress_image30](https://user-images.githubusercontent.com/45905048/69831085-20698700-1228-11ea-8549-7efad1ab0dc1.png)
