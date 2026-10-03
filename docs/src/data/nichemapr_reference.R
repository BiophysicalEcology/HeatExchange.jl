# Reference output of NicheMapR for the tutorials "An endotherm: metabolic rate" and "A human of many parts".
#
# Run once with NicheMapR installed; the docs read the CSV files this writes.
#   Rscript nichemapr_reference.R

library(NicheMapR)

# ── endoR_devel: the default 65 kg animal, no thermoregulation ──────────────
#
# DELTAR = 100 makes the exhaled air leave at lung temperature, as HeatExchange.jl does.
# DELTAR = 0 is the NicheMapR default (exhaled at air temperature).

TAs <- seq(0, 40, 2)
endo <- function(TA, DELTAR) {
  out <- endoR_devel(TA = TA, THERMOREG = 0, DELTAR = DELTAR)
  data.frame(air_temperature_C = TA, deltar = DELTAR,
             metabolic_W = out$enbal[, "QGEN"], evaporation_W = out$enbal[, "QEVAP"],
             convection_W = out$enbal[, "QCONV"], longwave_in_W = out$enbal[, "QIRIN"],
             longwave_out_W = out$enbal[, "QIROUT"],
             skin_C = out$treg[, "TSKIN_D"], fur_C = out$treg[, "TFA_D"], lung_C = out$treg[, "TLUNG"],
             fur_conductivity = out$treg[, "K_FUR_D"],
             respiratory_water_g_h = out$masbal[, "H2OResp_g"], cutaneous_water_g_h = out$masbal[, "H2OCut_g"],
             area_m2 = out$morph[, "AREA"], skin_area_m2 = out$morph[, "AREA_SKIN"])
}
endo_table <- do.call(rbind, c(lapply(TAs, endo, DELTAR = 100), lapply(TAs, endo, DELTAR = 0)))
write.csv(endo_table, "endoR_reference.csv", row.names = FALSE)

# ── HomoTherm: the default human in the cold, where no thermoregulation is invoked ──

parts <- c("head", "trunk", "arm", "leg")
homo <- function(TA) {
  out <- HomoTherm(TA = TA, VEL = 0.1, RH = 50)
  part_rows <- do.call(rbind, lapply(parts, function(p) {
    treg <- out[[paste0(p, ".treg")]]
    enbal <- out[[paste0(p, ".enbal")]]
    morph <- out[[paste0(p, ".morph")]]
    data.frame(air_temperature_C = TA, part = p,
               core_C = treg["T_CORE"], skin_dorsal_C = treg["TSKIN_D"], skin_ventral_C = treg["TSKIN_V"],
               surface_dorsal_C = treg["TCLO_D"], surface_ventral_C = treg["TCLO_V"],
               flesh_conductivity = treg["K_FLESH"], skin_wetness_pct = treg["PCTWET"],
               heat_generated_W = enbal["QMETAB"], evaporation_W = enbal["QEVAP"], convection_W = enbal["QCONV"],
               area_m2 = morph["AREA"], join_area_m2 = morph["AREA_JOIN"])
  }))
  b <- out$balance
  whole <- data.frame(air_temperature_C = TA, metabolic_W = b["QMETAB"], core_C = b["T_CORE"], lung_C = b["T_LUNG"],
                      skin_C = b["T_SKIN"], surface_C = b["T_CLO"], flesh_conductivity = b["K_FLESH"],
                      skin_wetness_pct = b["PCTWET"], respiration_evaporation_W = -b["QEVAP_RESP"],
                      respiration_convection_W = b["QCONV_RESP"], cutaneous_evaporation_W = -b["QEVAP_CUT"],
                      convection_W = -b["QCONV"], area_m2 = b["AREA"])
  list(parts = part_rows, whole = whole)
}
homo_out <- lapply(seq(-10, 30, 2), homo)
write.csv(do.call(rbind, lapply(homo_out, `[[`, "parts")), "homotherm_parts.csv", row.names = FALSE)
write.csv(do.call(rbind, lapply(homo_out, `[[`, "whole")), "homotherm_whole.csv", row.names = FALSE)
