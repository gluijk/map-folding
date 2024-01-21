# Folded maps
# www.overfitting.net
# https://www.overfitting.net/2024/01/plegando-papel-con-r.html


library(tiff)  # save 16-bit TIFF's
library(png)  # save 8-bit PNG's


hillshademap=function(DEM, dx=25, dlight=c(0, 2, 3), gamma=1) {
    # hillshademap() inputs DEM data and outputs a hillshade matrix
    #
    # DEM: digital elevation map matrix
    # dx: DEM resolution (cell size in the same units as elevation values)
    # dlight: lighting direction (3D vector defined from observer to light source):
    #   dlight=c(0, 2, 3)  # sunrise
    #   dlight=c(0, 0, 1)  # midday
    #   dlight=c(0,-2, 3)  # sunset
    # gamma: optional output gamma lift
    
    DIMY=nrow(DEM)
    DIMX=ncol(DEM)
    # If array turn its first dimension into a matrix
    if (!is.matrix(DEM)) DEM=matrix(DEM[,,1], nrow=DIMY, ncol=DIMX)
    
    dlightM=sum(dlight^2)^0.5
    
    # Vectorial product to calculate n (orthogonal vector)
    nx = 2*dx*(DEM[1:(DIMY-2), 2:(DIMX-1)] - DEM[3:DIMY,     2:(DIMX-1)])
    ny = 2*dx*(DEM[2:(DIMY-1), 1:(DIMX-2)] - DEM[2:(DIMY-1), 3:DIMX])
    nz = 4*dx^2
    nM = (nx^2 + ny^2 + nz^2)^0.5
    
    # Dot product to calculate cos(theta)
    dn = dlight[1]*nx + dlight[2]*ny + dlight[3]*nz  # (DIMY-2)x(DIMX-2) matrix
    
    # Reflectance (=cos(theta))
    hillshadepre=dn/(dlightM*nM)
    hillshadepre[hillshadepre<0]=0  # clip negative values
    
    # Add 1-pix 'lost' borders
    hillshademap=matrix(0, nrow=DIMY, ncol=DIMX)
    hillshademap[2:(DIMY-1), 2:(DIMX-1)]=hillshadepre
    rm(hillshadepre)
    hillshademap[c(1,DIMY),]=hillshademap[c(2,DIMY-1),]
    hillshademap[,c(1,DIMX)]=hillshademap[,c(2,DIMX-1)]
    
    return(hillshademap^(1/gamma))
}


###########################################################

map=readPNG("michiganhorizons.png")

DIMY=nrow(map)
DIMX=ncol(map)


# Folding parameters
NFOLDSY=3  # number of vertical folds
NFOLDSX=1  # number of horizontal folds
ENHANCE=2  # folding function peakness (higher=sharp, lower=progressive)
DELTA=floor(min(DIMY/(NFOLDSY+1), DIMX/(NFOLDSX+1))/2)

# Folding function (LENGTH = DELTA + 1 + DELTA)
LENGTH=2*DELTA+1  # total length of folding function
x=seq(from=0, to=pi/2, length.out=DELTA+1)  # half folding function
f=((1-cos(x)))^ENHANCE
f=c(f, rev(f[1:DELTA]))  # mirror folding function
plot(seq(from=-1, to=1, length.out=LENGTH), f,
     type='l', col='red', main='Folding function (1 - cos(x))^4',
     xlab='', ylab='', ylim=c(0,1), xaxp=c(-1, 1, 2), yaxp=c(0, 1, 1))
abline(h=0, v=0, lty='dotted')


# Create folding "DEM"'s
# They are arrays ranging [-1,1] with 0 as neutral value (elevation)
demy=array(0, c(DIMY, DIMX))  # blank fold sheet
HALF=round(DIMX/2)
for (i in 1:NFOLDSY) {
    rowfold=round(DIMY/(NFOLDSY+1))*i
    demy[(rowfold-DELTA):(rowfold+DELTA), 1:HALF       ]=f*(-1)^(i+1)
    demy[(rowfold-DELTA):(rowfold+DELTA), (HALF+1):DIMX]=f*(-1)^i
}

demx=t(demy)*0
for (i in 1:NFOLDSX) {
    colfold=round(DIMX/(NFOLDSX+1))*i
    demx[(colfold-DELTA):(colfold+DELTA), 1:DIMY]=f*(-1)^i
}
demx=t(demx)

# Refine demy to soften encounter between vertical and horizontal folds
for (i in 1:DIMY) {
    demy[i, (colfold-DELTA):(colfold+DELTA)]=
        demy[i, (colfold-DELTA):(colfold+DELTA)]*(1-f)
}

# Save normalizedd ([-1,1] -> [0,1]) folding "DEM"'s
writePNG((demy+1)/2, "demy.png")
writePNG((demx+1)/2, "demx.png")


# Calculate grayscale hillshades
dx=1/40  # "DEM"'s range [-1,1], 1/40 is a good tradeoff
# lower dx=more contrast but can cause clipping in the hillshade

# hillshade y
MIX=0.85  # two light sources are mixed to fill dark areas a bit
hilly=hillshademap(demy, dx=dx, dlight=c(-1, -2, 3))
hillfilly=hillshademap(demy, dx=dx, dlight=c(-1, -3, 2))
hilly=hilly*MIX+hillfilly*(1-MIX)

# hillshade x
MIX=0.85  # two light sources are mixed to fill dark areas a bit
hillx=hillshademap(demx, dx=dx, dlight=c(-1, -2, 3))
hillfillx=hillshademap(demx, dx=dx, dlight=c(-1, -3, 2))
hillx=hillx*MIX+hillfillx*(1-MIX)

# Check hillshade ranges (avoid clipping)
print(paste0("hilly (min/med/max): ", round(min(hilly),3), ' / ',
             round(median(hilly),3), ' / ', round(max(hilly),3)))
print(paste0("hillx (min/med/max): ", round(min(hillx),3), ' / ',
             round(median(hillx),3), ' / ', round(max(hillx),3)))


# Save and display hillshades

# hillshade y
writeTIFF(hilly^2, "hillshadey.tif",
          bits.per.sample=16, compression="LZW")
# Display hillshade
image(t(hilly[nrow(hilly):1,]), useRaster=TRUE,
      col=c(gray.colors(256, start=0, end=1, gamma=0.5)),
      asp=nrow(hilly)/ncol(hilly), axes=FALSE)

# hillshade x
writeTIFF(hillx^2, "hillshadex.tif",
          bits.per.sample=16, compression="LZW")
# Display hillshade
image(t(hillx[nrow(hillx):1,]), useRaster=TRUE,
      col=c(gray.colors(256, start=0, end=1, gamma=0.5)),
      asp=nrow(hillx)/ncol(hillx), axes=FALSE)


# Apply folding to map (multiply map * hillshade) and save
# Both hillshades are averaged and then multiplied by the map
med=median(hilly)  # median(hilly)=median(hillx)
gamma=log(0.5*med)/log(0.5)  # this gamma restores 0.5*med into 0.5
writeTIFF( ( map * replicate(3, hilly+hillx)/2) ^(1/gamma),
           "foldedmap.tif", bits.per.sample=16, compression="LZW")
