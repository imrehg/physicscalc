"""
All the wires I need to use and then some.
"""


# class Wire(object):
#     diameter = 0
#     resistance = 0
#     name

#     def __init__(self, dia, res, name):
#         """
#         Initialize with wire parameters

#         dia: diameter (m)
#         res: resistance (Ohm / m)
#         """
#         self.diameter = dia
#         self.resistance = res
#         self.name = name

#     def getDia(self):
#         return self.diameter

#     def getResistance(self):
#         return self.resistance

#     def __str__(self):
#         print self.name

# # Pre-create a bunch of American Wire Gauge wires
# AWG1 = Wire(7.348e-3, 0.4066e-3, "AWG1")
# AWG2 = Wire(6.544e-3, 0.5127e-3, "AWG2")
# AWG3 = Wire(5.827e-3, 0.6465e-3, "AWG3")
# AWG4 = Wire(5.189e-3, 0.8152e-3, "AWG4")
# AWG5 = Wire(4.621e-3, 1.028e-3, "AWG5")
# AWG6 = Wire(4.115e-3, 1.296e-3, "AWG6")
# AWG7 = Wire(3.665e-3, 1.634e-3, "AWG7")
# AWG8 = Wire(3.264e-3, 2.061e-3, "AWG8")
# AWG9 = Wire(2.906e-3, 2.599e-3, "AWG9")
# AWG10 = Wire(2.588e-3, 3.277e-3, "AWG10")

# Pre-create a bunch of American Wire Gauge wires
AWG1 = (7.348e-3, 0.4066e-3, "AWG1")
AWG2 = (6.544e-3, 0.5127e-3, "AWG2")
AWG3 = (5.827e-3, 0.6465e-3, "AWG3")
AWG4 = (5.189e-3, 0.8152e-3, "AWG4")
AWG5 = (4.621e-3, 1.028e-3, "AWG5")
AWG6 = (4.115e-3, 1.296e-3, "AWG6")
AWG7 = (3.665e-3, 1.634e-3, "AWG7")
AWG8 = (3.264e-3, 2.061e-3, "AWG8")
AWG9 = (2.906e-3, 2.599e-3, "AWG9")
AWG10 = (2.588e-3, 3.277e-3, "AWG10")

# Wire collection
AWG = [AWG1, AWG2, AWG3, AWG4, AWG5, AWG6, AWG7, AWG8, AWG9, AWG10]
