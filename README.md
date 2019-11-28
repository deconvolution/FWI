# FWI
Acoustic Full-waveform inversion (FWI) for isotropic media, where the wave equation is solved in the frequency domain with perfectly mathced layers (PML) and absorbing boundary condition (ABC) by finite difference method. It takes about one hour for each iteration as the wave equation is solved in the frequency domain, which is much slower than the time domain solution. The progress of each step will be generated in the command window, as below

![WeChat Screenshot_20191128214948](https://user-images.githubusercontent.com/45905048/69831281-08463780-1229-11ea-89de-6fdb4e72bc65.png)
\
Forward and adjoint wavefields are cross-correlated to calculate the gradient for velocity correction.

# Generate true recordings
Run 'generate_true_rec.m' to generate true recordings 'true_rect.mat' based on a toy model 'trueModel'. Below is our true velocity model (stolen from the courst TPG4155 at NTNU).
![true velocity](https://user-images.githubusercontent.com/45905048/68906553-f81c5b80-0744-11ea-90ab-1384d10d7f28.jpg)
# Run FWI with true recordings
Run 'FWI.mat' to process FWI based on the true recordings to reconstruct the velocity model. The starting model is given below.
![starting velocity](https://user-images.githubusercontent.com/45905048/68906711-855fb000-0745-11ea-812b-d576b0eb66ec.jpg)
Below is the final result after 30 times of interations. The final velocity model is close to our true velocity model. The salt structure on the left of the true model is not well reconstructed. The reason is that the high velocity of the salt bends the P-wave considerably, making it hard for P-wave penetrate the inner side. 
![progress_image30](https://user-images.githubusercontent.com/45905048/69831085-20698700-1228-11ea-8549-7efad1ab0dc1.png)
# Appendix: FWI progress
## Iteration 1
![progress_image1](https://user-images.githubusercontent.com/45905048/69831122-4bec7180-1228-11ea-96ea-f8493f29ddf2.png)
## Iteration 5
![progress_image5](https://user-images.githubusercontent.com/45905048/69831143-5c9ce780-1228-11ea-94d0-dcf456292c58.png)
## Iteration 10
![progress_image10](https://user-images.githubusercontent.com/45905048/69831155-69b9d680-1228-11ea-845b-2fe8027502cf.png)
## Iteration 15
![progress_image15](https://user-images.githubusercontent.com/45905048/69831162-73dbd500-1228-11ea-9aa2-8a4d9fa53bc2.png)
## Iteration 20
![progress_image20](https://user-images.githubusercontent.com/45905048/69831171-7d653d00-1228-11ea-82b8-7a560596ae7e.png)
## Iteration 25
![progress_image25](https://user-images.githubusercontent.com/45905048/69831180-89e99580-1228-11ea-908c-0aded1c15606.png)
## Iteration 30
![progress_image30](https://user-images.githubusercontent.com/45905048/69831186-9241d080-1228-11ea-90fe-0e9112de720c.png)
