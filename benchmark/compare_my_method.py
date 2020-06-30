import os
import sys
# sys.path.insert(0, 'common')
import benchmark.single.scripts.plotroutines as case1
# import benchmark.field.scripts.plotroutines as case4
# import single.scripts.plotroutines as case1
# import field.scripts.plotroutines as case4
import matplotlib.pyplot as plt
import numpy as np

def calculate_match(x, y, ls, lowerpercentile, upperpercentile):
    x_eval = np.linspace(max(min(x), min(ls)), min(max(x), max(ls)), 100)
    y_eval = np.interp(x_eval, x, y)
    lower_eval = np.interp(x_eval, ls, lowerpercentile)
    upper_eval = np.interp(x_eval, ls, upperpercentile)
    inbetween = np.logical_and(np.less_equal(lower_eval, y_eval), np.greater_equal(upper_eval, y_eval))
    return sum(inbetween.astype(int))

def evaluate(x, y, case, plot_id, ref_index=0, line_or_fracture_id=0):
    ax = plt.gca()

    if case == 1:

        places_and_methods = {
            "UiB": ["TPFA", "MPFA", "MVEM", "RT0"],
            "USTUTT": ["MPFA", "TPFA\_Circ"],
            "LANL": ["MFD"],
            "NCU\_TW": ["Hybrid\_FEM"],
            "UNICE\_UNIGE": ["VAG\_Cont", "HFV\_Cont", "VAG\_Disc", "HFV\_Disc"],
            "ETHZ\_USI": ["FEM\_LM"],
            "UNICAMP": ["Hybrid\_Hdiv"],
            "UNIL\_USI": ["FE\_AMR\_AFC"],
            "INM": ["EDFM"],
            "DTU": ["FEM\_COMSOL"]
        };
        # os.chdir("single/scripts")
        os.chdir("../../benchmark/single/scripts")

        if plot_id == "a":
            (ls, lp, up) = case1.plot_percentiles(str(ref_index), case1.id_p_matrix,
                                                  places_and_methods, ax)
        elif plot_id == "b":
            (ls, lp, up) = case1.plot_percentiles(str(ref_index), case1.id_c_matrix,
                                                  places_and_methods, ax)
        elif plot_id == "c":
            (ls, lp, up) = case1.plot_percentiles(str(ref_index), case1.id_c_fracture,
                                                  places_and_methods, ax)
        else:
            print("Error. Invalid plot id " + plot_id + " provided.")
            os.chdir("../..")
            sys.exit(1)

    # elif case == 4:

    #     places_and_methods = {
    #         "UiB": ["TPFA", "MPFA", "MVEM", "RT0"],
    #         "USTUTT": ["MPFA", "TPFA\_Circ"],
    #         "LANL": ["MFD"],
    #         "UNICE\_UNIGE": ["VAG\_Cont", "HFV\_Cont", "VAG\_Disc", "HFV\_Disc"],
    #         "ETHZ\_USI": ["FEM\_LM"],
    #         "UNICAMP": ["Hybrid\_Hdiv"],
    #         "DTU": ["FEM\_COMSOL"]
    #     };
    #     os.chdir("field/scripts")

    #     if plot_id == "a":
    #         if line_or_fracture_id == 2:
    #             line_or_fracture_id = 0
    #         (ls, lp, up) = case4.plot_percentiles(str(line_or_fracture_id), places_and_methods, ax)
    #     else:
    #         print("Error. Invalid plot id " + plot_id + " provided.")
    #         os.chdir("../..")
    #         sys.exit(1)

    else:

        print("Error. Invalid case id " + str(case) + " provided.")
        sys.exit(1)

    match = calculate_match(x, y, ls, lp, up)

    ax.plot(x, y, 'C1', linewidth=3)
    ax.set_title("match: " + str(match))
    # ax.legend(["your method", "published spread"])
    ax.legend(["This method", "Berre et al (2020)"])

    # os.chdir("../..")

    return plt