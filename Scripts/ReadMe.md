# Scripts

## Analayses

### NMDS ordination
- saved in `NMDS`
- Generates a 3D NMDS ordination for invertebrate communnity data
  - Using bray-curtis dissimilarity matrix
  - 3D solution reached after 85 iterations with a stress of 0.2042451
  - Scaling: centring, PC rotation, halfchange scaling 

### Cluster analysis and perMANOVA
- saved in `Cluster` 
- Trying to do an iterative approach to agglomerative hierarchical clustering
  - use ISA to find appropriate number of groups
  - should use flexible beta (-0.25) linkages
  - did it in PCORD to get it done quickly (l o l )
 - perMANVOA
    - used `adonis` (Anderson 2001 approach)
    - 999 permutations
    - using Bray-Curtis dissimilarity matrix
 - used `betadisper` to look at homogeneity among groups - actually pretty good!
