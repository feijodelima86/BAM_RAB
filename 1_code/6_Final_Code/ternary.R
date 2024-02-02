library("Ternary")

# If using your own data, set
# abc <- [Three-column matrix containing a, b, c coordinates of points]

# We'll use random data for this example.

abc <- rbind(a, b, c)

response <- d

# Now we must start a plot, to define the coordinate system
par(mar = rep(0.2, 4))
TernaryPlot(alab = "Cd_TRD", blab = "TURNOVER", clab = "OM")

# Convert measured points to XY
xy <- TernaryToXY(abc)

# Use an inverse distance weighting to interpolate between measured points
Predict <- function(predXY) {
  Distance <- function(a, b) {
    apply(a, 2, function(pt) sqrt(colSums((pt - b) ^ 2)))
  }
  dists <- Distance(xy, predXY)
  id <- 1 / dists
  idw <- id / rowSums(id)
  
  # Return:
  colSums(response * t(idw))
}

# Predict at triangle centres
tri <- TriangleCentres(resolution = 15L)

# Adjust the resolution to suit your own dataset

# Now we interpolate between our known values to generate a colour for each
# of our tiles

predicted <- Predict(tri[1:2, ])
map <- rbind(x = tri["x", ], y = tri["y", ], z = predicted,
             down = tri["triDown", ])

# Place a semitransparent colour fill over grid lines:
ColourTernary(map)

# Calculate contours
PredictABC <- function(a, b, c) Predict(TernaryToXY(rbind(a, b, c)))
TernaryContour(PredictABC, resolution = 100L, legend = 6, bty = "n")

# Mark the points at which we took measurements
TernaryPoints(abc, pch = 20, col = "black")

