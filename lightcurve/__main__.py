import numpy as np
import dearpygui.dearpygui as dpg

from __init__ import Orbit

DATA = np.arange(0, 10, 0.001)


class GUI:
    def __init__(self):
        self.orbit = Orbit(3 / 215, 0.1, 0 * np.pi, 0, 0)
        self.orbit.init_primaries(1, 1, 2, 1, 1, 10)
        self.orbit.init_distance(121)
        self.curve = [self.orbit.app_mag(t) for t in DATA]
        dpg.create_context()
        self.init_window()
        dpg.create_viewport(title="LightCurve", width=800, height=800)
        dpg.setup_dearpygui()
        dpg.show_viewport()
        dpg.start_dearpygui()
        dpg.destroy_context()

    def init_window(self):
        with dpg.window(label="Plot", tag="win"):
            dpg.add_slider_float(
                label="a",
                callback=lambda _, a: [
                    self.orbit.__setattr__("a", a),
                    self.__setattr__("curve", [self.orbit.app_mag(t) for t in DATA]),
                    dpg.set_value("series_tag", [DATA, self.curve]),
                ],
                min_value=0.001,
                max_value=10,
                default_value=self.orbit.a,
            )
            dpg.add_slider_float(
                label="e",
                callback=lambda _, e: [
                    self.orbit.__setattr__("e", e),
                    self.__setattr__("curve", [self.orbit.app_mag(t) for t in DATA]),
                    dpg.set_value("series_tag", [DATA, self.curve]),
                ],
                min_value=0,
                max_value=1,
                default_value=self.orbit.e,
                clamped=True,
            )
            dpg.add_slider_float(
                label="i",
                callback=lambda _, i: [
                    self.orbit.__setattr__("i", i),
                    self.orbit.init_conversion(),
                    self.__setattr__("curve", [self.orbit.app_mag(t) for t in DATA]),
                    dpg.set_value("series_tag", [DATA, self.curve]),
                ],
                min_value=0,
                max_value=2 * np.pi,
                default_value=self.orbit.i,
                clamped=True,
            )
            dpg.add_slider_float(
                label="omega",
                callback=lambda _, omega: [
                    self.orbit.__setattr__("omega", omega),
                    self.orbit.init_conversion(),
                    self.__setattr__("curve", [self.orbit.app_mag(t) for t in DATA]),
                    dpg.set_value("series_tag", [DATA, self.curve]),
                ],
                min_value=0,
                max_value=2 * np.pi,
                default_value=self.orbit.omega,
                clamped=True,
            )
            dpg.add_slider_float(
                label="D",
                callback=lambda _, D: [
                    self.orbit.__setattr__("D", D),
                    self.__setattr__("curve", [self.orbit.app_mag(t) for t in DATA]),
                    dpg.set_value("series_tag", [DATA, self.curve]),
                ],
                min_value=0.001,
                max_value=1000,
                default_value=self.orbit.D,
                clamped=True,
            )
            dpg.add_slider_float(
                label="R1",
                callback=lambda _, R1: [
                    self.orbit.__setattr__("R1", R1),
                    self.__setattr__("curve", [self.orbit.app_mag(t) for t in DATA]),
                    dpg.set_value("series_tag", [DATA, self.curve]),
                ],
                min_value=0.001,
                max_value=10,
                default_value=self.orbit.R1,
            )
            dpg.add_slider_float(
                label="R2",
                callback=lambda _, R2: [
                    self.orbit.__setattr__("R2", R2),
                    self.__setattr__("curve", [self.orbit.app_mag(t) for t in DATA]),
                    dpg.set_value("series_tag", [DATA, self.curve]),
                ],
                min_value=0.001,
                max_value=10,
                default_value=self.orbit.R2,
            )
            dpg.add_slider_float(
                label="L1",
                callback=lambda _, L1: [
                    self.orbit.__setattr__("L1", L1),
                    self.__setattr__("curve", [self.orbit.app_mag(t) for t in DATA]),
                    dpg.set_value("series_tag", [DATA, self.curve]),
                ],
                min_value=0.001,
                max_value=10,
                default_value=self.orbit.L1,
            )
            dpg.add_slider_float(
                label="L2",
                callback=lambda _, L2: [
                    self.orbit.__setattr__("L2", L2),
                    self.__setattr__("curve", [self.orbit.app_mag(t) for t in DATA]),
                    dpg.set_value("series_tag", [DATA, self.curve]),
                ],
                min_value=0.001,
                max_value=10,
                default_value=self.orbit.L2,
            )

            # create plot
            with dpg.plot(label="Lightcurve", height=500, width=750):
                dpg.add_plot_legend()

                dpg.add_plot_axis(dpg.mvXAxis, label="Time (years)", tag="x_axis")
                dpg.add_plot_axis(
                    dpg.mvYAxis, label="Magnitude", tag="y_axis", invert=True
                )

                dpg.add_line_series(DATA, self.curve, parent="y_axis", tag="series_tag")
                dpg.set_axis_limits_constraints("x_axis", 0, 3)


if __name__ == "__main__":
    GUI()
