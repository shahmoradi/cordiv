# This R script uses the package Rpdb to extract the CA coordinates of a given pdb file which will be then used to generate a two dimensional Voronoi diagram.
# Last updated by Amir Shhamoradi, 10:15 AM, Thursday Feb 19 2015, WilkeLab, ICMB, UT Austin

# USES library Rpdb for reading parsing pdb file
# USES library deldir for generating 2D Voronoi diagram

getcrd = function(file,atom,crd)
{
  pdb <- read.pdb(file)
  crd = pdb$atoms[pdb$atoms$elename == atom, c('x1','x2','x3')]
  return(crd)
}

cacrd = getcrd(file='../../structures/RepairPDB_1LBA_A.pdb',atom='CA')

vt = deldir(cacrd$x1,cacrd$x2)  # generate Voronoi tessellationa dn Delaunay trangulation
filename = '../figures/voronoi_diagram.pdf'
pdf( filename, width=4.5, height=4, useDingbats=FALSE )
par( mai=c(0.65, 0.65, 0.1, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 )
plot( vt
    , wlines='tess'
    , pch=19:1
    , cex=0.3
    , lty = 1
    , col = c(34,1)
    , xlab = 'X coordinate [ Angstroms ]'
    , ylab = 'Y coordinate [ Angstroms ]'
    , showrect = TRUE
    )
graphics.off()

