library("dplyr")
library("tidyr")
library("devtools")

locs <- read.table(system.file("extdata","Stationinfo.dat",package = "STRbook"),
                   col.names = c("id","lat","lon"))

Times <- read.table(system.file("extdata","Times_1990.dat",package = "STRbook"),
                    col.names = c("julian","year","month","day"))

Tmax <- read.table(system.file("extdata","Tmax_1990.dat", package = "STRbook"))

#rename the Tmax cols to the same names as locs cols
names(Tmax) <- locs$id

#combine the Times and Tmax df
Tmax <- cbind(Times, Tmax)
head(names(Tmax),10)

#use gather to put the dataframe into long format. This function takes the data as first argument, 
#the key-value pair, and then the next arguments are the names of any columns to exclude as values.

Tmax_long <- gather(Tmax, id, z, -julian, -year, -month, -day)
head(Tmax_long)

Tmax_long$id <- as.integer(Tmax_long$id)

#filter out data with values less than -9998
nrow(Tmax_long)

Tmax_long <- filter(Tmax_long, !(z <= -9998))
nrow(Tmax_long)

#use mutate to add a col to df
Tmax_long <- mutate(Tmax_long, proc = "Tmax")
head(Tmax_long)

#load the data directly from STRbook
data(Tmin_long, package = "STRbook")
data(TDP_long, package = "STRbook")
data(Precip_long, package = "STRbook")

#construct final dataframe in long format by simply concatenating all these(rowwise) together
#suing the function rbind
NOAA_df_1990 <- rbind(Tmax_long, Tmin_long, TDP_long, Precip_long)

#the function group_by creates a grouped dataframe, while summarise does an operation on each group within 
#grouped dataframe
#find the mean value for each variable in each year
summ <- group_by(NOAA_df_1990, year, proc) %>%
  summarise(mean_proc = mean(z))


#find the number of days on which it did not rain at each station in June of every year
NOAA_precip <- filter(NOAA_df_1990, proc == "Precip" & month == 6)
summ <- group_by(NOAA_precip, year, id) %>%
  summarise(days_no_precip = sum(z == 0))
head(summ)

median(summ$days_no_precip)

NOAA_df_sorted <- arrange(NOAA_df_1990, julian, id)

head(NOAA_df_sorted)

#merge locs and NOAA_df_1990 using the function left_join
NOAA_df_1990 <- left_join(NOAA_df_1990, locs, by = "id")

#constructs a space-wide dataframe of max temp, with each row denoting a different date and each col containing
#data z from a specific station id
Tmax_long_sel <- select(Tmax_long, julian, id, z)
Tmax_wide <- spread(Tmax_long_sel, id, z)
dim(Tmax_wide)

M <- select(Tmax_wide, -julian) %>% as.matrix()

#Working with Spatiotemporal data classes
library("sp")
library("spacetime")

#construct date in year-month-day format using function paste
# The with( ) function applys an expression to a datase
NOAA_df_1990$date <- with(NOAA_df_1990,
                          paste(year, month, day, sep = "-"))
head(NOAA_df_1990$date, 4)

#convert date from character to a Date object
NOAA_df_1990$date <- as.Date(NOAA_df_1990$date)
class(NOAA_df_1990$date)

#stConstruct class requires spatial and temporal coordinates
Tmax_long2 <- filter(NOAA_df_1990, proc == "Tmax")
STObj <- stConstruct(x = Tmax_long2,                    #data set
                     space = c("lon","lat"),            #spatial fields
                     time = "date")                     #time fields
class(STObj)

#use function STIDF, require to specify the spatial part as an object of class Spatial from package sp
spat_part <- SpatialPoints(coords = Tmax_long2[, c("lon","lat")])
temp_part <- Tmax_long2$date
STObj2 <- STIDF(sp = spat_part,
                time = temp_part,
                data = select(Tmax_long2, -date, -lon, -lat))
class(STObj2)

#constructing a STFDF object
spat_part <- SpatialPoints(coords = locs[, c("lon","lat")])
temp_part <- with(Times, 
                  paste(year, month, day, sep = "-"))
temp_part <- as.Date(temp_part)

Tmax_long3 <- gather(Tmax, id, z, -julian, -year, -month, -day)

Tmax_long3$id <- as.integer(Tmax_long3$id)
Tmax_long3 <- arrange(Tmax_long3, julian, id)

#confirming that the spatial ordering is the correct one 
all(unique(Tmax_long3$id) == locs$id)
#construct STFDF
STObj3 <- STFDF(sp = spat_part,
                time = temp_part,
                data = Tmax_long3)

class(STObj3)

#equip STObj3 with a coordinate reference system
proj4string(STObj3) <- CRS("+proj=longlat + ellps=WGS84")

#replace the missing values with NAs
STObj3$z[STObj3$z == -9999] <- NA

#Visualization
library("animation")
library("dplyr")
library("ggplot2")
library("gstat")
library("maps")
library("STRbook")

set.seed(1)
#load data set and take a subset using function filter
data("NOAA_df_1990", package = "STRbook")
Tmax <- filter(NOAA_df_1990,             #subset the data
               proc == "Tmax" &          #only max temperature
                 month %in% 5:9 &        #May to Spetember
                 year == 1993)           #year of 1993

Tmax %>% select(lon, lat, date, julian, z) %>% head()

#create a variable t that is equal to 1 when julian == 728050 and increases by 1 for each day in the record

Tmax$t <- Tmax$julian - 728049       #create a new time variable

Tmax %>% select(lon, lat, date, julian, z, t) %>% head()

#Spatial Plots
Tmax_1 <- subset(Tmax, t %in% c(1, 15, 30))          #extract data
NOAA_plot <- ggplot(Tmax_1) +         #plot points
  geom_point(aes(x = lon, y = lat,    #lon and lat
                 colour = z),         #attribute color
             size = 2) +              #make all points larger
  col_scale(name = "degF") +          #attach color scale
  xlab("Longitude (deg)") +           #x-axis label
  ylab("Latitude (deg)") +          #y-axis label
  geom_path(data = map_data("state"),
            aes(x = long, y = lat, group = group)) +
  facet_grid(~date) +                 #facet by time
  coord_fixed(xlim = c(-100, -75),
              ylim = c(25, 50)) + 
  theme_bw()
print(NOAA_plot)

#Plot point-referenced data
#use geom_tile for regular lattice data
#use geom_polygon for irregular lattice data

#load the BEA data
data("BEA", package = "STRbook")
head(BEA %>% select(-Description), 3)

#The data contains boundary points for the counties
data("MOcounties", package = "STRbook")
head(MOcounties %>% select(long, lat, NAME10), 3)
#to plot the boundary of the first county one can simply type:
County1 <- filter(MOcounties, NAME10 == "Clark, MO")
plot(County1$long, County1$lat)

#use left join to add geo data to income data
MOcounties <- left_join(MOcounties, BEA, by = "NAME10")

#use geom_polyon to display BEA income data as spatial polyon
#use geom_path to draw the county boundaries
#use group argument to identify which points correspond to which county
g1 <- ggplot(MOcounties) +
  geom_polygon(aes(x = long, y = lat,      #county boundary
                  group = NAME10,         #county group
                  fill = log(X1970))) +   #log of income
  geom_path(aes(x = long, y = lat,        #county boundary
                group = NAME10)) +        #county group 
  fill_scale(limits = c(7.5, 10.2), 
             name = "log($)") +
  coord_fixed() + ggtitle("1970") +       #annotations
  xlab("x (m)") + ylab("y (m)") + theme_bw()
print(g1)    

#Time Series Plots
#choose 10 at random and extract the data associated with these 10 stations from the data set. 
UIDs <- unique(Tmax$id)                      #extract IDs
UIDs_sub <- sample(UIDs, 10)                 #sample 10 IDs
Tmax_sub <- filter(Tmax, id %in% UIDs_sub)   #subset data

TmaxTS <- ggplot(Tmax_sub) +
  geom_line(aes(x = t, y = z)) +         #line plot of z against t
  facet_wrap(~id, ncol = 5) +            #facet by station
  xlab("Day number  (days)") +           #x label
  ylab("Tmax (degF)") +                  #y label
  theme_bw() +                           #BW theme
  theme(panel.spacing = unit(1, "lines")) # facet spacing
print(TmaxTS)  

#Hovmoller Plots
#a 2D space-time visualization where space is collapsed(projected or averaged) onto one dimension; the second dimension denotes time.

lim_lat <- range(Tmax$lat)          #latitude range
lim_t <- range(Tmax$t)
lat_axis <- seq(lim_lat[1],         #latitude axis
                lim_lat[2],
                length = 25)
t_axis <- seq(lim_t[1], 
              lim_t[2],
              length = 100)
#generate a regular grid using the function expand.grid
lat_t_grid <- expand.grid(lat = lat_axis,
                          t = t_axis)
#next, we need to associate each station's latitudinal coordinate with the closest one on the grid.
#this can be done by finding the distance from the station's latitudinal coordinate to each point of the grid,
#finding which gridpoint is the closest

Tmax_grid <- Tmax
dist <- abs(outer(Tmax$lat, lat_axis, "-"))
Tmax_grid$lat <- lat_axis[apply(dists, 1, which.min)]


  
  
  
  
  







































