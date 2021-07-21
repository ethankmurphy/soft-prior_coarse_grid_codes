# soft-prior_coarse_grid_codes

The project has 2 example meshes to help show the functionality of general soft-prior coarse grid, coarse-to-fine mapping, and Laplace smoothing regularization construction codes. The codes rely on a mesh and functions based on surface triangulations representing subdomains to work. There is some dependency on codes in S:\digihisto\Ethan\EKM_utility

Run construct_inv_grids_test.m to see how the soft-prior code works.

The main function is constr_coarsesub_grid_genfuncs.m, it makes the coarse grids that include the soft-prior information, but how the functions are called in construct_inv_grids_test.m are important to get it to work correctly. 

There are additionally two functions that help 1) merge_grids_maps_Lmats.m and 2) check_softprior_cgrid_c2f_L.m, which do what their names imply

The meshes are in the subfolder dat
* sphere_1inc_r10_2_h0o8_0o3_0o6: Sphere with 1 inclusion
* sphere_shell_1inc_r10_8_2_h0o8_0o3_0o6_0o6: sphere with a large spherical inner subdomain and an inclusion within this subdomain

Other files used:
*simpplot: From distmesh (unchanged)
* plot_msh_elecs_only: Plots only the electrodes
* ndrm_laplacian_regmatrix: in NDRM
* eval_signed_dist: in Ethan\EKM_utility
* get_tcs: in Ethan\EKM_utility
* vis_reg_connections: in Ethan\EKM_utility
