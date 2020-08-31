# PhD Chapter 3
## Invertebrate communities


Comparing the invertebrate communities among invaded (Phragmites), uninvaded (cattail, meadow) and herbicide-treated (restore) sites in Long Point. There were 9 sites in each vegetation type and they spanned a water depth gradient. Were collected in two years, 2017 and 2018, over the field season (4 collections in 2017, 6 in 2018). Will combine them since differences between years are not of interest. 

### Analyses

1. NMDS ordination

    - To visualize differences among the vegetation types; likely will keep both years together
    - Data was relativized by column maximum, and removed rare Families (anything with 2 or fewer occurrences)
    - Performed ordination using metaMDS with bray-curtis distances (in script `Explore` right now)
        - extracted the coordinates using scores()
        - used envfit to find the Families that had an r2 > 0.2 for axis 1,2 and axis 1,3 and pulled them out
