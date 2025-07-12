library(sf)
library(tidyverse)
library(ggplot2)
library(countrycode)
library(igraph)
library(geosphere)

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
      warning("Stop")
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


###################### Whole europe


path <- paste("./Data_task_45/EGM_2019_SHP_20190312/DATA/FullEurope/RailrdC.shp", sep ="")
# read the GIS file
shape.file <- as.data.frame(st_read(path, quiet = T))
euNodes <- processDataFrame(shape.file)

path <- "./Data_task_45/EGM_2019_SHP_20190312/DATA/FullEurope/RailrdL.shp"
df <- as.data.frame(st_read(path, quiet = TRUE))

# filter on operational only lines
df |> filter(EXS == 28) -> df
# Now extract the coordinates only
coordsLink <- as.data.frame(st_coordinates(df$geometry))
coordsLink <- coordsLink |> 
  group_by(L1) |> 
  summarise(
    X_start = first(X),
    Y_start = first(Y),
    X_end   = last(X),
    Y_end   = last(Y)
  ) |> 
  ungroup()
# apply the findNearestStation
res <- lapply(1:nrow(coordsLink), function(i) findNearestStation(coordsLink[i,], euNodes, 70000))
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
  findUpstreamChain(downstreams[d,], nonConnected, 70000) -> chain
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
finalDF_clean -> edges
edges <- edges %>% distinct(to, from, .keep_all = TRUE)
edges |> filter(from != to) -> edges

pdf("EU.pdf")
g <- graph_from_data_frame(
  d = edges[,c("to", "from")],
  vertices = euNodes,
  directed = FALSE
)
layout_coords <- cbind(euNodes$longitude, euNodes$latitude)

plot(g,
     layout = layout_coords,
     vertex.size = 0.5,
     vertex.label = NA,
     vertex.color = "red",
     vertex.frame.color = NA,       # Nessun contorno visibile
     edge.arrow.size = 0.05,
     edge.width = 0.5,
     main = "EU")

dev.off()


# And write everything to memory
dir.create("./data/EU", recursive = TRUE, showWarnings = FALSE)

# write vertices
write.table(euNodes,"./data/EU/vertices.txt", row.names = F)
write.table(edges[,c("from", "to")],"./data/EU/edges.txt", row.names = F)



############################ Degree distribution ############################
# Compute the degree distribution
degree_distribution <- degree(g, mode = "all")
# Create a data frame for plotting
degree_df <- data.frame(degree = degree_distribution)
# Plot the degree distribution
pdf("DDEU.pdf")
p<-ggplot(degree_df, aes(x = degree)) +
  geom_histogram(binwidth = 1, fill = "blue", color = "black", alpha = 0.7) +
  labs(title = paste("Degree Distribution for", country.code),
       x = "Degree",
       y = "Frequency") +
  theme_minimal() 
print(p)
dev.off()