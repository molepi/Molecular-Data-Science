# install.packages("clusterSim") #Install package
# library("clusterSim") #load package
# source('c:/DBindex.r')
# 
# Example
# output = function(x, md, hc, min_nc, max_nc)
# output$clusters
# output$res
# 

DBindex <- function(x, md, hc, min_nc, max_nc){

# make array to store DB-index
res = array(0, c(max_nc-min_nc+1, 2))

# make number of clusters to compute DB-index
res[,1] = min_nc:max_nc  

# Start with NULL
clusters = NULL 

# RUN for all pre-defined clusters
for (nc in min_nc:max_nc)
{
	cl2 = cutree(hc, k=nc) # Take nc clusters
	res[nc-min_nc+1,2] = DB = index.DB(x,cl2,d=md,centrotypes="medoids")$DB # Compute DB
	clusters = rbind(clusters, cl2) # Store results
}


# Show optimal number of clusters
print(paste("minimum Davies-Bouldin index for",(min_nc:max_nc)[which.min(res[,2])],"clusters=",min(res[,2])))
# print(clusters[which.min(res[,2]),])

# Make plot
plot(res[,2], type="o", pch=0, xlab="Number of clusters", ylab="DB", xaxt="n")
axis(1, c(min_nc:max_nc))


# Define output
list(clusters=clusters,res=res)

}