'''
2016 Amy Campbell
Example script using staNMF class. Reproduces Wu et al. 2016
staNMF for drosophila spatial expression data by
running staNMF for # principal patterns K between 15 and 30
'''
print("1")
import staNMF as st
print("2")
# Construct NMF object with custom output folder name
# Note: filename = 'example' uses the included
# 'updatedexpressionpatterns.csv'  and parameter 'folderID'
# is not required
exampleNMF = st.staNMF(filename='example',
                       folderID="DrosophilaExpression")
print("3")
# Run NMF and output factorization solutions over the range
# of K's
exampleNMF.runNMF()
print("4")
# Calculate instability for each K
exampleNMF.instability(k1=15, k2=30)

# Plot instability against K with a custom title
exampleNMF.plot(dataset_title="Example staNMF in Drosophila"
                "Spatial Expression Data 8-13")
