library(NicheMapR)

ellipsoid_input = list(
posture = 4.5,
mass = 0.5,
density = 1000,
coreT = 37,
furdepth = 5,
furcond = 0.04,
O2eff = 0.2,
stress = 0.6,
windspd = 0.1,
rh = 50,
Q10 = 3,
basal = NA,
basmult = 1)
airT <- seq(5, 45)

endo <- ellipsoid(posture = ellipsoid_input$posture,
                mass = ellipsoid_input$mass,
                density = ellipsoid_input$density,
                coreT = ellipsoid_input$coreT,
                furdepth = ellipsoid_input$furdepth,
                furcond = ellipsoid_input$furcond,
                O2eff = ellipsoid_input$O2eff,
                stress = ellipsoid_input$stress,
                airT = airT,
                windspd = ellipsoid_input$windspd,
                rh = ellipsoid_input$rh,
                Q10 = ellipsoid_input$Q10,
                basal = ellipsoid_input$basal,
                basmult = ellipsoid_input$basmult)
ellipsoid_output <- as.data.frame(endo)

ellipsoid_input$basal <- 0.0

write.csv(unlist(ellipsoid_input), file = '../data/ellipsoid_input.csv')
write.csv(airT, file = '../data/ellipsoid_air_temperature.csv')
write.csv(ellipsoid_output, file = '../data/ellipsoid_output.csv')
