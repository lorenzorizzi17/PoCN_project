library(sf)
library(tidyverse)
library(ggplot2)
library(countrycode)
library(igraph)
library(geosphere)

# Browse through the filesystem and see which folders are available
countryCodesISO2 <- system("ls ./Data_task_45/EGM_2019_SHP_20190312/DATA/Countries/", intern = TRUE)
cat("Found the following country codes: \n")
cat(countryCodesISO2)
cat("\n")

# Select only valid countries where railroad data is available
validCountries <- c()
for(cc in countryCodesISO2){
  system(paste("ls ./Data_task_45/EGM_2019_SHP_20190312/DATA/Countries/", cc,"/Rail* 2>/dev/null", sep = ""), intern = TRUE) -> output
  if (length(output)==0){
    cat("\nHaven't found Rail data for", countrycode(cc, origin = "iso2c", destination = "country.name"))
  } else {
    validCountries <- c(validCountries, cc)
  }
}

cat("\n\nValid countries where data has been found:\n")
cat(validCountries)


## This function reads and extracts the dataset from GIS file
extractStationDataFrame <- function(country.code){
  # absolute path
  path <- paste("./Data_task_45/EGM_2019_SHP_20190312/DATA/Countries/", country.code, "/RailrdC.shp", sep ="")
  # read the GIS file
  shape.file <- as.data.frame(st_read(path, quiet = T))
  return(shape.file)
}

# Process the vertices dataframe to get it into the right format
processDataFrame <- function(shape.file){
  coords <- st_coordinates(shape.file$geometry)
  shape.file <- cbind(shape.file, coords)

  nodeID <- integer(length(shape.file$NAMA1))
  nodeID[order(shape.file$NAMN1)] = 1:length(shape.file$NAMN1)
  df <- data.frame(
    nodeID = nodeID,
    nodeLabel = shape.file$NAMN1,
    latitude = coords[,2],
    longitude = coords[,1],
    country_name = countrycode(shape.file$ICC, origin = "iso2c", destination = "country.name"),
    country_ISO3 = countrycode(shape.file$ICC, origin = "iso2c", destination = "iso3c")
  )
  return(df)
}

extractRailRoadsDataFrame <- function(country.code){
  path <- paste0("./Data_task_45/EGM_2019_SHP_20190312/DATA/Countries/", country.code, "/RailrdL.shp")
  shape.file <- st_read(path, quiet = TRUE)
  return(shape.file)
}


#####################################################################
## given a link, it finds the stations (of origin and destination) ##
## closest to its endpoints. If none is found, returns NA          ##
#####################################################################
findNearestStation <- function(link, dfCountry, minimumDistanceStation){
  origin <- as.numeric(link[c("X_start", "Y_start")])
  dest <- as.numeric(link[c("X_end", "Y_end")])
  
  origin_dists <- distHaversine(cbind(dfCountry$longitude, dfCountry$latitude),
                                matrix(origin, nrow = 1))
  dest_dists <- distHaversine(cbind(dfCountry$longitude, dfCountry$latitude),
                              matrix(dest, nrow = 1))
  minOrigin <- min(origin_dists)
  minDestin <- min(dest_dists)
  origin_label <- if(minOrigin < minimumDistanceStation) dfCountry$nodeID[which.min(origin_dists)] else NA
  dest_label <- if(minDestin < minimumDistanceStation) dfCountry$nodeID[which.min(dest_dists)] else NA
  
  return(list(origin = origin_label, destination = dest_label))
}

#####################################################################
### Given a downstream node, tries to build the chain and         ###
### connect its root node (i.e. upstreams), eventually parsing    ###
### intermediates                                                 ###
#####################################################################
findUpstreamChain <- function(downstream, potentialUpstreams, minimumDistanceRails) {
  currentLink <- downstream
  chain <- list()

  repeat {  # repeat until un upstream node is reached
    # coords of current node
    currentCoords <- matrix(c(currentLink$X_start, currentLink$Y_start), nrow=1)
    upstreamCoords <- as.matrix(potentialUpstreams[, c("X_end", "Y_end")])
    
    if (
      nrow(potentialUpstreams) == 0 ||
      any(is.na(currentCoords)) || ncol(currentCoords) != 2 || nrow(currentCoords) == 0 ||
      any(is.na(upstreamCoords)) || ncol(upstreamCoords) != 2 || nrow(upstreamCoords) == 0
    ) {
      warning("Coordinate non valide o upstream vuoto, interrompo la ricerca")
      break
    }


    distances <- distHaversine(upstreamCoords, currentCoords)
    idx_within_threshold <- which(distances < minimumDistanceRails)

    # Haven't found any connection. ISsue a warning and discard the downstream node
    if (length(idx_within_threshold) == 0) {
      warning("Not found a connection")
      break
    }
    
    nearest_idx <- idx_within_threshold[which.min(distances[idx_within_threshold])]
    # Get the candidate upstream node
    upstream <- potentialUpstreams[nearest_idx[1], , drop = FALSE]

    # Build the dataframe row to eventually add
    new_row <- data.frame(
      L1 = downstream$L1,
      X_start = upstream$X_start,
      Y_start = upstream$Y_start,
      X_end = currentLink$X_end,
      Y_end = currentLink$Y_end,
      from = upstream$from,
      to = currentLink$to
    )
    
    # add to the end of the list the new row
    chain[[length(chain) + 1]] <- new_row

    # If the upstream node is a station, we can stop here. Break from repeat
    if (!is.na(upstream$from)) {
      break
    }
    # Otherwise, repeat considering the upstream node as the currentLink

    currentLink <- upstream
    # and delete the already computed
    potentialUpstreams <- potentialUpstreams[-nearest_idx[1], , drop = FALSE]
  } #end repeat

  if (length(chain) > 0) {
    return(do.call(rbind, chain))
  } else {
    return(NULL)
  }
}

######################################
## Put all those functions together ##
######################################
produceEdgeDataframe <- function(country.code, minStation, minRails){
  # Extract one of those dataframes
  df <- as.data.frame(extractRailRoadsDataFrame(country.code))
  # filter on operational only lines
  df |> filter(EXS == 28) -> df
  # Now extract the coordinates only
  coordsLink <- as.data.frame(st_coordinates(df$geometry))
  coordsLink <- coordsLink |> group_by(L1) |>
    summarise(
      X_start = first(X),
      Y_start = first(Y),
      X_end = last(X),
      Y_end = last(Y)
    ) |> ungroup()
  # apply the findNearestStation
  res <- lapply(1:nrow(coordsLink), function(i) findNearestStation(coordsLink[i,], dataFramesStations[[country.code]], minStation))
  coordsLink$from <- sapply(res, `[[`, "origin")
  coordsLink$to <- sapply(res, `[[`, "destination")
  
  # reject lines connecting one station to itself
  coordsLink <- coordsLink |> filter(!( !is.na(to) & !is.na(from) & to == from ))

  # Edges that are definitive (they already connect two stations)
  finalDF <- coordsLink |> filter( !is.na(to) & !is.na(from) )
  # Or downstream edges
  downstreams <- coordsLink |> filter(!is.na(to) & is.na(from))
  # Or non connected edges (potential upstreams)
  nonConnected <- coordsLink |> filter(is.na(to))

  # iterator over the downstreams
  for(d in 1:nrow(downstreams)){
    # get the chain
    findUpstreamChain(downstreams[d,], nonConnected, minRails) -> chain
    # and if succesfull, add it to the definitive df
    if (!is.null((chain))){
      toAdd <- data.frame(L1 = chain[1,"L1"],
                        X_start = chain[ nrow(chain), "X_start"],
                        Y_start = chain[ nrow(chain), "Y_start"],
                        X_end = chain[ 1, "X_start"],
                        Y_end = chain[1, "Y_start"],
                        from = chain[ nrow(chain), "from"],
                        to = chain[1, "to"])
      finalDF <- rbind(finalDF, toAdd)
      }
    }
  
  finalDF_clean <- finalDF[!is.na(finalDF$from) & !is.na(finalDF$to), ]
}

# We are finally read to run the algorithm

# A list of dataframes, one for each country
dataFramesStations <- list()
for (cc in validCountries){
  df <- extractStationDataFrame(cc)
  df <- processDataFrame(df)
  dataFramesStations[[cc]] <- df 
}


# Choose the country
country.code = validCountries[11]
print(paste("\n Selected: ", country.code))

# Extract the vertices dataframes and the edge dataframes ...
vertices = dataFramesStations[[country.code]]
edges = produceEdgeDataframe(country.code, minStation = 50000, minRails = 50000)
# remove double links
edges <- edges %>% distinct(to, from, .keep_all = TRUE)

# Draw the graph
pdf(paste(country.code, ".pdf"))

g <- graph_from_data_frame(
  d = edges[,c("to", "from")],
  vertices = vertices,
  directed = FALSE
)
layout_coords <- cbind(vertices$longitude, vertices$latitude)

plot(g,
     layout = layout_coords,
     vertex.size = 1.5,
     vertex.label = NA,
     vertex.color = "red",
     vertex.frame.color = NA,       
     edge.arrow.size = 0.05,
     edge.width = 0.5,
     main = country.code)
dev.off()
print("\n Graph built \n")

## Now the actual drawing ...

pdf(paste("Real",country.code ,".pdf"))

extractStationDataFrame(country.code) -> nodes_sf
extractRailRoadsDataFrame(country.code) -> edges_sf
nodes_sf <- st_as_sf(nodes_sf)
edges_sf <- st_as_sf(edges_sf)

ggplot() +
  geom_sf(data = edges_sf, color = "gray50", size = 0.3) +       # le ferrovie
  geom_sf(data = nodes_sf, color = "red", size = 1.5) +          # le stazioni
  theme_minimal() +
  labs(title = country.code) +
  coord_sf()

dev.off()

# And write everything to memory
dir.create(paste("./data/", country.code, sep = ""), recursive = TRUE, showWarnings = FALSE)

# write vertices
write.table(vertices,paste("./data/", country.code, "/vertices.txt", sep = ""), row.names = F)
write.table(edges[,c("from", "to")],paste("./data/", country.code, "/edges.txt", sep = ""), row.names = F)
