"""Passive, responsive charts which do no rendering while their tab is hidden."""

from nicegui.elements.plotly import Plotly


class VisiblePlotly(Plotly, component="visible_plotly.js"):
    """Keep NiceGUI's figure serialization; resize only visible chart panels.

    The inherited Plotly import registers its JS dependency. This component is
    for the flight deck's passive figures; no Python-side Plotly events are used.
    """

