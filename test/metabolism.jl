using HeatExchange
using Unitful

mass = 1.0u"kg"
rq = 0.8

Q_metab = metabolic_rate(Kleiber(), mass)
Q_metab = metabolic_rate(McKechnieWolf(), mass)

O2_consumption = u"ml/hr"(Joules_to_O2(Typical(), Q_metab, rq))
O2_consumption = u"ml/hr"(Joules_to_O2(Q_metab))

u"W"(O2_to_Joules(Typical(), O2_consumption, rq))
u"J/hr"(O2_to_Joules(1.0u"ml/hr"))

u"ml/hr"(Joules_to_O2(Kleiber1961(), Q_metab, rq))
u"W"(O2_to_Joules(Kleiber1961(), O2_consumption, rq))