#Moran's I and Geary's C


library(spatstat)                              #Spatial statistical analysis by many authors
library(dplyr)                                 #A fast, consistent tool for working with data frame like objects by Hadley Wickam et al.
library(smerc)                                 #Implements statistical methods for analyzing the counts of areal data, with a 
                                               #focus on detection of spatial clusters and clustering. The package emphasis on spatial scan methods!
library(sp)                                    #Classes and methods for spatial data
library(ape)                                   #Package for analyses of phylogenetics and evolution
library(sf)                                    #Work with vector data in R
library(maptools)                              #Conversion from sf to owin requires the use of the maptools package

library(tmap)

# Set the working directory, I always like to do this so I don't lose files and to simplify subsequent read and writes
setwd("C:/Users/yx2957/Documents/Research/Induced Seismicity/R Codes")                           # choose your local working directory

# Read the data table from a comma delimited file - data on GitHub/GeostatsGuy/GeoDataSets
mydata = read.csv("Result_Pts.csv")        # read in comma delimited data file with all injection and all events

# Let's visualize the first several rows of our data so we can make sure we successfully loaded it
head(mydata) # show the first several rows of a data table in the console
mydata

# Check out the summary statistics for each column
summary(mydata)                                # summary statistics for the multivariate data file

# Subset the dataframe
mydata_event = subset(mydata, TYPE == 0)
mydata_inj = subset(mydata, TYPE == 1)

crds_event <- cbind(x=mydata_event$SURFX, y=mydata_event$SURFY)
P1 <- Polygon(crds_event)

Pls <- SpatialPolygons(list(P1))



# Assign the spatial extents 
xlim = c(1733000,2250000); ylim = c(152300,718900)

# Plots the points and calculate and plot the convex hull (clock wise)
plot(mydata_event$SURFX, mydata_event$SURFY, cex = 0.5,xlab = 'X (ft)',ylab='Y(ft)',xlim = xlim, ylim = ylim, main='Events and Convex Hull',asp=1.0)
#switched the shape from convex to rectangle
hpts_event = ripras(mydata_event$SURFX, mydata_event$SURFY, shape="convex")
plot(hpts_event,add=TRUE)



#creates an object that representing a point pattern dataset in the 2D plane, with a window
pp_event <- ppp(mydata_event$SURFX,mydata_event$SURFY,window = hpts_event, unitname=c("meters","meters"))
# Notice that we have just pased the locations (x,y) and the extents of the box aligned with the x and y coordinates and units


# Let's check for and remove duplicates
pp_unique_event = unique(pp_event,warn=TRUE)
pp_unique_event
summary(pp_unique_event)
plot(pp_unique_event, cex = 0.5, pch = "+")

#convert ppp object to a spatialpoints object, we can use the following from the maptools package:
pp_unique_event_sp <- as(pp_unique_event, "SpatialPoints")
class(pp_unique_event_sp)
plot(pp_unique_event_sp)



#use quadrats to divides window into rectangle quadrats and count the number in each quadrats
#The return value is a table which printed neatly. The return value is also a member of the special class "quadratcount". Plotting the object
#will display the quadrats, annotated by their counts. 
rec_quadrats <- quadrats(pp_unique_event, 10, 10)
quadrats_count <- quadratcount(pp_unique_event, 10, 10)
plot(rec_quadrats, main = "Events quadrat plot")
plot(quadrats_count, add=T)

noisy <- im(quadrats_count)
plot(noisy)

########Compute Moran's I###########
b <- quantile(quadrats_count, probs = (0:3)/3)
slope_break <- cut(quadrats_count, breaks = 5, labels = 1:5)
#A tessellation is a collection of disjoint spatial regions that fit together o form a larger spatial region. 
p.owin <- as.im(slope_break)
v <- tess(image = slope_break)

#Kolmogorov-Smirnov test of CSR. test for the complete spatial randomness
# I am using density of events for cdf.test. NEED TO CHANGE TO PORE PRESSURE AS COVARIATE!!!
ks <- cdf.test(pp_unique_event, ds, "ks")
plot(ks)
##########################################################################################################

#Nearest neighbor distance (G function): the distance of each point to its nearest neighbor
#Empty space distance (F function): the distance from a fixed reference location to the nearest data point
#Pairwise distance (K function): the distance between all pairs of points

#Random labeling: given the locations X, marks are conditionally independent and identically ditributed(i.i.d)
#Independence of components: the sub-process Xm for each mark m are independent point processes
#CSR and independence: the locations are a uniform Poisson point process and the marks are i.i.d


#Convert the datafrmae to a spatial object
#p.sf <- st_as_sf(mydata, coords = c("SURFX","SURFY"))

#convert to an sf object; owin can be coerced to an sf object via the st_as_sf
pp.sf <- st_as_sf(pp_unique_event_sp)
#st_as_sf(s.sp)

#Converting an sf object to a spatial object/ SpatialpolygonsDataFrame
pp.sp <- as(pp.sf, 'Spatial')

par(mai=c(0,0,0,0))
plot(pp.sp)
xy <- coordinates(pp.sp)

library(spdep)
w <- poly2nb(pp.sp, row.names=pp.sp$Id)



#Exporting to different data file formats
st_write(pp.sf, "mydata_out.shp", driver = "ESRI Shapefile", delete_layer=TRUE)  # create to a shapefile

plot(pp.sf)






































































