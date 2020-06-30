linestyle = {"UiB" : {"MVEM": "-", "TPFA": "-", "MPFA": "-", "RT0": ":"},
             "USTUTT" : {"MPFA": "-", "TPFA\_Circ": "-", "reference": "-"},
             "UNICE\_UNIGE" : {"HFV\_Cont": "-", "HFV\_Disc": ":", "VAG\_Cont": "-", "VAG\_Disc": ":"},
             "NCU\_TW" : {"Hybrid\_FEM": "--"},
             "LANL" : {"MFD\_Tet": "--", "MFD\_Hex": "--", "MFD\_Ex": "--", "MFD": "--"},
             "ETHZ\_USI" : {"FEM\_LM": "--"},
             "UNIL\_USI" : {"FE\_AMR\_AFC": "--"},
             "UNICAMP" : {"Hybrid\_Hdiv": "--"},
             "INM" : {"EDFM": "--"},
             "DTU" : {"FEM\_COMSOL": "--"},
            }

color = {"UiB" : {"MVEM": "C0", "TPFA": "C1", "MPFA": "C2", "RT0": "C3"},
         "USTUTT" : {"MPFA": "C4", "TPFA\_Circ": "C5", "reference": "black"},
         "UNICE\_UNIGE" : {"HFV\_Cont": "C6", "HFV\_Disc": "C7", "VAG\_Cont": "C8", "VAG\_Disc": "C9"},
         "NCU\_TW" : {"Hybrid\_FEM": "C0"},
         "LANL" : {"MFD\_Tet": "C1", "MFD\_Hex": "C2", "MFD\_Ex": "C1", "MFD": "C2"},
         "ETHZ\_USI" : {"FEM\_LM": "C3"},
         "UNIL\_USI" : {"FE\_AMR\_AFC": "C4"},
         "UNICAMP" : {"Hybrid\_Hdiv": "C5"},
         "INM" : {"EDFM": "C6"},
         "DTU" : {"FEM\_COMSOL": "C7"}
        }

# Returns the label used for time
def getTimeLabel(unit='s'):
    label = "$t \, [\mathrm{"
    label += unit
    label += "}]$"
    return label

# Returns the label used for arc lengths
def getArcLengthLabel():
    return "$\mathrm{arc} \, \mathrm{length} \, [m]$"

# Returns the label used for piezometric head
def getHeadLabel(dimension):
    # label = "$h_"
    label = "$p_"
    label += str(dimension)
    label += " \, [\mathrm{m}]$"
    return label

# Returns the label used for concentrations
def getConcentrationLabel(dimension):
    label = "$c_"
    label += str(dimension)
    label += " \, [\mathrm{m}^{-3}]$"
    return label

# Returns the label used for avered concentrations
def getAveragedConcentrationLabel(dimension):
    label = "$\overline{c_"
    label += str(dimension)
    label += "} \, [\mathrm{m}^{-3}]$"
    return label
