# Supporting Code and Data to replicate RRmorph case studies

**Reference paper**: Melchionna, M., Castiglione, S., Girardi, G., Serio, C., Esposito, A., Mondanaro, A., Profico, A., Sansalone, G., & Raia, P.(2024).RRmorph—a new R package to map phenotypic evolutionary rates and patterns on 3D meshes. Communications Biology 7: 1009.doi:10.1038/s42003-024-06710-8

**Repository content:**

**CASE STUDIES - DATA FILES**

**case study 1 ratemap.rda**

The file contains:

1. data: data.frame with information on species and group
2. mshapeS: 3D surface of the consensus shape
3. pairedLM: landmarks and semilandmarks indices
4. refmat: landmarks and semilandmarks configurations used for the interpolation
5. refsur: 3D meshes used for the interpolation
6. slid: not-superimposed landmarks and semilandmarks configurations
7. tree: phylogenetic tree

**case study 2 convmap.rda**

The file contains:

1. data: data.frame with information on species and group
2. mshapeS: 3D surface of the consensus shape
3. pairedLM: landmarks and semilandmarks indices
4. refmat: landmarks and semilandmarks configurations to be used for the interpolation
5. refsur: 3D meshes to be used for the interpolation
6. slid: not-superimposed landmarks and semilandmarks configurations
7. tree: phylogenetic tree

**case study 3 convmap skull and endocast.rda**

The file contains:

1. data: data.frame with information on species and group
2. mshapeS_e: 3D surface of the endocasts consensus shape
3. mshapeS_c: 3D surface of the skulls consensus shape
4. pairedLM_e: landmarks and semilandmarks indices of endocasts configuration
5. pairedLM_c: landmarks and semilandmarks indices of skulls configuration
6. refmat: landmarks and semilandmarks configurations to be used for the interpolation
7. refsur: 3D meshes to be used for the interpolation
8. slid_e: not-superimposed landmarks and semilandmarks endocast configurations
9. slid_c: not-superimposed landmarks and semilandmarks skull configurations
10. tree: phylogenetic tree

### SUPPORTING R CODE

The scripts to run all data as in the case studies detailed above, and to plot figures
