using HeatExchange
using Unitful

Q_metab = 1.0u"W"
O2_consumption = 1.0u"ml/hr"
rq = 0.8

u"ml/hr"(Joules_to_O2(Typical(), Q_metab))
u"ml/hr"(Joules_to_O2(Q_metab))
u"J/hr"(O2_to_Joules(Typical(), O2_consumption))
u"J/hr"(O2_to_Joules(1.0u"ml/hr"))

u"ml/hr"(Joules_to_O2(Kleiber1961(), Q_metab, rq))
u"J/hr"(O2_to_Joules(Kleiber1961(), O2_consumption, rq))
