library(NicheMapR)

# plate
GEOMETRY <- 0
shape_b <- 3
shape_c <- 2/3
AMASS <- 0.04
ORIENT <- 0 # TODO
SHP <- c(1, shape_b, shape_c)
CUSTOMGEOM <- c(10.4713, 0.688, 0.425, 0.85, 3.798, 0.683, 0.694, 0.743) # TODO
ANDENS <- 1000
SKINW <- 0.001 # TODO
SKINT <- 0 # TODO
RINSUL <- 0 # TODO
PTCOND <- 0.1 # TODO
PMOUTH <- 0.05 # TODO
PANT <- 1 # TODO
plate_input <- c(GEOMETRY, shape_b, shape_c, AMASS, ANDENS, SKINW, SKINT, RINSUL, PTCOND, PMOUTH, PANT)
GEOM_out <- GEOM_ecto(AMASS = AMASS,
                      GEOMETRY = GEOMETRY,
                      ORIENT = ORIENT,
                      SHP = c(1, shape_b, shape_c),
                      CUSTOMGEOM = CUSTOMGEOM,
                      ANDENS = ANDENS,
                      SKINW = SKINW,
                      SKINT = SKINT,
                      RINSUL = RINSUL,
                      PTCOND = PTCOND,
                      PMOUTH = PMOUTH,
                      PANT = PANT)
write.csv(plate_input, '../data/plate_in.csv')
write.csv(unlist(GEOM_out), '../data/plate_out.csv')

# cylinder
GEOMETRY <- 1
shape_b <- 3
shape_c <- 2/3 # redundant
AMASS <- 0.04
ORIENT <- 0 # TODO
SHP <- c(1, shape_b, shape_c)
CUSTOMGEOM <- c(10.4713, 0.688, 0.425, 0.85, 3.798, 0.683, 0.694, 0.743) # TODO
ANDENS <- 1000
SKINW <- 0.001 # TODO
SKINT <- 0 # TODO
RINSUL <- 0 # TODO
PTCOND <- 0.1 # TODO
PMOUTH <- 0.05 # TODO
PANT <- 1 # TODO
cylinder_input <- c(GEOMETRY, shape_b, shape_c, AMASS, ANDENS, SKINW, SKINT, RINSUL, PTCOND, PMOUTH, PANT)
GEOM_out <- GEOM_ecto(AMASS = AMASS,
                      GEOMETRY = GEOMETRY,
                      ORIENT = ORIENT,
                      SHP = c(1, shape_b, shape_c),
                      CUSTOMGEOM = CUSTOMGEOM,
                      ANDENS = ANDENS,
                      SKINW = SKINW,
                      SKINT = SKINT,
                      RINSUL = RINSUL,
                      PTCOND = PTCOND,
                      PMOUTH = PMOUTH,
                      PANT = PANT)
write.csv(cylinder_input, '../data/cylinder_in.csv')
write.csv(unlist(GEOM_out), '../data/cylinder_out.csv')


# ellipsoid
GEOMETRY <- 2
shape_b <- 3
shape_c <- 2/3
AMASS <- 0.04
ORIENT <- 0 # TODO
SHP <- c(1, shape_b, shape_c)
CUSTOMGEOM <- c(10.4713, 0.688, 0.425, 0.85, 3.798, 0.683, 0.694, 0.743) # TODO
ANDENS <- 1000
SKINW <- 0.001 # TODO
SKINT <- 0 # TODO
RINSUL <- 0 # TODO
PTCOND <- 0.1 # TODO
PMOUTH <- 0.05 # TODO
PANT <- 1 # TODO
ellipsoid_input <- c(GEOMETRY, shape_b, shape_c, AMASS, ANDENS, SKINW, SKINT, RINSUL, PTCOND, PMOUTH, PANT)
GEOM_out <- GEOM_ecto(AMASS = AMASS,
                      GEOMETRY = GEOMETRY,
                      ORIENT = ORIENT,
                      SHP = c(1, shape_b, shape_c),
                      CUSTOMGEOM = CUSTOMGEOM,
                      ANDENS = ANDENS,
                      SKINW = SKINW,
                      SKINT = SKINT,
                      RINSUL = RINSUL,
                      PTCOND = PTCOND,
                      PMOUTH = PMOUTH,
                      PANT = PANT)
write.csv(ellipsoid_input, '../data/ellipsoid_in.csv')
write.csv(unlist(GEOM_out), '../data/ellipsoid_out.csv')

# desert iguana
GEOMETRY <- 3
shape_b <- 3 # redundant
shape_c <- 2/3 # redundant
AMASS <- 0.04
ORIENT <- 0 # TODO
SHP <- c(1, shape_b, shape_c)
CUSTOMGEOM <- c(10.4713, 0.688, 0.425, 0.85, 3.798, 0.683, 0.694, 0.743) # TODO
ANDENS <- 1000
SKINW <- 0.001 # TODO
SKINT <- 0 # TODO
RINSUL <- 0 # TODO
PTCOND <- 0.1 # TODO
PMOUTH <- 0.05 # TODO
PANT <- 1 # TODO
iguana_input <- c(GEOMETRY, shape_b, shape_c, AMASS, ANDENS, SKINW, SKINT, RINSUL, PTCOND, PMOUTH, PANT)
GEOM_out <- GEOM_ecto(AMASS = AMASS,
                      GEOMETRY = GEOMETRY,
                      ORIENT = ORIENT,
                      SHP = c(1, shape_b, shape_c),
                      CUSTOMGEOM = CUSTOMGEOM,
                      ANDENS = ANDENS,
                      SKINW = SKINW,
                      SKINT = SKINT,
                      RINSUL = RINSUL,
                      PTCOND = PTCOND,
                      PMOUTH = PMOUTH,
                      PANT = PANT)
write.csv(iguana_input, '../data/iguana_in.csv')
write.csv(unlist(GEOM_out), '../data/iguana_out.csv')

# leopard frog
GEOMETRY <- 4
shape_b <- 3 # redundant
shape_c <- 2/3 # redundant
AMASS <- 0.04
ORIENT <- 0 # TODO
SHP <- c(1, shape_b, shape_c)
CUSTOMGEOM <- c(10.4713, 0.688, 0.425, 0.85, 3.798, 0.683, 0.694, 0.743) # TODO
ANDENS <- 1000
SKINW <- 0.001 # TODO
SKINT <- 0 # TODO
RINSUL <- 0 # TODO
PTCOND <- 0.1 # TODO
PMOUTH <- 0.05 # TODO
PANT <- 1 # TODO
frog_input <- c(GEOMETRY, shape_b, shape_c, AMASS, ANDENS, SKINW, SKINT, RINSUL, PTCOND, PMOUTH, PANT)
GEOM_out <- GEOM_ecto(AMASS = AMASS,
                      GEOMETRY = GEOMETRY,
                      ORIENT = ORIENT,
                      SHP = c(1, shape_b, shape_c),
                      CUSTOMGEOM = CUSTOMGEOM,
                      ANDENS = ANDENS,
                      SKINW = SKINW,
                      SKINT = SKINT,
                      RINSUL = RINSUL,
                      PTCOND = PTCOND,
                      PMOUTH = PMOUTH,
                      PANT = PANT)
write.csv(frog_input, '../data/frog_in.csv')
write.csv(unlist(GEOM_out), '../data/frog_out.csv')
