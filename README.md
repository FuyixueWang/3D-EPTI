# 3D-EPTI
Example Matlab codes and data for 3D-EPTI reconstruction using low-rank subspace method.

Related works:
 1) "Ultrafast high-resolution multi-parametric brain MRI using 3D Echo Planar Time-resolved Imaging",Fuyixue Wang and Zijing Dong et al., 2021.
 2) "Echo planar time-resolved imaging (EPTI)", Fuyixue Wang et al., MRM, 2019.
 3) "EPTI with Subspace Reconstruction and Optimized Spatiotemporal Encoding", Zijing Dong et al., MRM, 2020.
 4) "T2 shuffling: Sharp, multicontrast, volumetric fast spin-echo imaging", Jonathan Tamir et al., MRM, 2017.

The example raw k-space data used for testing (script_3DEPTI_subspaceRecon_example_CPU.m) can be download at:
https://doi.org/10.6084/m9.figshare.14558154

The data for basis generation is in 'acq_params.mat'.
'Generated_bases.mat' is the generated bases.

This example presents a CPU implementation of the reconstruction based on BART toolbox: https://mrirecon.github.io/bart

Another example Matlab implementation without using BART can be found at https://github.com/zijingd/VFA-EPTI

Further acceleration can be achieved by GPU.

Please contact us if you have any questions about our work.
Fuyixue Wang <fwang18@mgh.harvard.edu>, Zijing Dong <zijingd@mit.edu> Feb/2021
